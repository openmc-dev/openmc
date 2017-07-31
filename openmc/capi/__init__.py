"""
This module provides bindings to C functions defined by OpenMC shared library.
When the :mod:`openmc` package is imported, the OpenMC shared library is
automatically loaded. Calls to the OpenMC library can then be via functions or
objects in the :mod:`openmc.capi` subpackage, for example:

.. code-block:: python

    openmc.capi.init()
    openmc.capi.run()
    openmc.capi.finalize()

"""

from contextlib import contextmanager
from ctypes import CDLL, c_int, c_int32, c_double, POINTER
import sys
from warnings import warn

from numpy.ctypeslib import as_array
import pkg_resources


# Determine shared-library suffix
if sys.platform == 'darwin':
    _suffix = 'dylib'
else:
    _suffix = 'so'

# Open shared library
_filename = pkg_resources.resource_filename(
    __name__, '_libopenmc.{}'.format(_suffix))
try:
    _dll = CDLL(_filename)
    _available = True
except OSError:
    warn("OpenMC shared library is not available from the Python API. This "
         "means you will not be able to use openmc.capi to make in-memory "
         "calls to OpenMC.")
    _available = False


class GeometryError(Exception):
    pass


def _error_code(s):
    """Get error code corresponding to global constant."""
    return c_int.in_dll(_dll, s).value


def _error_handler(err, func, args):
    """Raise exception according to error code."""
    if err == _error_code('e_out_of_bounds'):
        raise IndexError('Array index out of bounds.')

    elif err == _error_code('e_cell_not_allocated'):
        raise MemoryError("Memory has not been allocated for cells.")

    elif err == _error_code('e_cell_invalid_id'):
        raise KeyError("No cell exists with ID={}.".format(args[0]))

    elif err == _error_code('e_cell_not_found'):
        raise GeometryError("Could not find cell at position ({}, {}, {})"
                            .format(*args[0]))

    elif err == _error_code('e_nuclide_not_allocated'):
        raise MemoryError("Memory has not been allocated for nuclides.")

    elif err == _error_code('e_nuclide_not_loaded'):
        raise KeyError("No nuclide named '{}' has been loaded.")

    elif err == _error_code('e_nuclide_not_in_library'):
        raise KeyError("Specified nuclide doesn't exist in the cross "
                       "section library.")

    elif err == _error_code('e_material_not_allocated'):
        raise MemoryError("Memory has not been allocated for materials.")

    elif err == _error_code('e_material_invalid_id'):
        raise KeyError("No material exists with ID={}.".format(args[0]))

    elif err == _error_code('e_tally_not_allocated'):
        raise MemoryError("Memory has not been allocated for tallies.")

    elif err == _error_code('e_tally_invalid_id'):
        raise KeyError("No tally exists with ID={}.".format(args[0]))

    elif err < 0:
        raise Exception("Unknown error encountered (code {}).".format(err))


def calculate_volumes():
    """Run stochastic volume calculation"""
    _dll.openmc_calculate_volumes()


def finalize():
    """Finalize simulation and free memory"""
    _dll.openmc_finalize()


def find(xyz, rtype='cell'):
    """Find the cell or material at a given point

    Parameters
    ----------
    xyz : iterable of float
        Cartesian coordinates of position
    rtype : {'cell', 'material'}
        Whether to return the cell or material ID

    Returns
    -------
    int or None
        ID of the cell or material. If 'material' is requested and no
        material exists at the given coordinate, None is returned.
    int
        If the cell at the given point is repeated in the geometry, this
        indicates which instance it is, i.e., 0 would be the first instance.

    """
    # Set second argument to openmc_find
    if rtype == 'cell':
        r_int = 1
    elif rtype == 'material':
        r_int = 2
    else:
        raise ValueError('Unknown return type: {}'.format(rtype))

    # Call openmc_find
    uid = c_int32()
    instance = c_int32()
    _dll.openmc_find((c_double*3)(*xyz), r_int, uid, instance)
    return (uid.value if uid != 0 else None), instance.value


def hard_reset():
    """Reset tallies, timers, and pseudo-random number generator state."""
    _dll.openmc_hard_reset()


def init(intracomm=None):
    """Initialize OpenMC

    Parameters
    ----------
    intracomm : mpi4py.MPI.Intracomm or None
        MPI intracommunicator

    """
    if intracomm is not None:
        # If an mpi4py communicator was passed, convert it to an integer to
        # be passed to openmc_init
        try:
            intracomm = intracomm.py2f()
        except AttributeError:
            pass
        _dll.openmc_init(c_int(intracomm))
    else:
        _dll.openmc_init(None)


def keff():
    """Return the calculated k-eigenvalue and its standard deviation.

    Returns
    -------
    tuple
        Mean k-eigenvalue and standard deviation of the mean

    """
    k = (c_double*2)()
    _dll.openmc_get_keff(k)
    return tuple(k)


def plot_geometry():
    """Plot geometry"""
    _dll.openmc_plot_geometry()


def reset():
    """Reset tallies and timers."""
    _dll.openmc_reset()


def run():
    """Run simulation"""
    _dll.openmc_run()


@contextmanager
def run_in_memory(intracomm=None):
    """Provides context manager for calling OpenMC shared library functions.

    This function is intended to be used in a 'with' statement and ensures that
    OpenMC is properly initialized/finalized. At the completion of the 'with'
    block, all memory that was allocated during the block is freed. For
    example::

        with openmc.capi.run_in_memory():
            for i in range(n_iters):
                openmc.capi.reset()
                do_stuff()
                openmc.capi.run()

    Parameters
    ----------
    intracomm : mpi4py.MPI.Intracomm or None
        MPI intracommunicator

    """
    init(intracomm)
    yield
    finalize()


# Set argument/return types
if _available:
    _dll.openmc_calculate_volumes.restype = None
    _dll.openmc_finalize.restype = None
    _dll.openmc_find.argtypes = [
        POINTER(c_double*3), c_int, POINTER(c_int32), POINTER(c_int32)]
    _dll.openmc_find.restype = c_int
    _dll.openmc_find.errcheck = _error_handler
    _dll.openmc_hard_reset.restype = None
    _dll.openmc_init.argtypes = [POINTER(c_int)]
    _dll.openmc_init.restype = None
    _dll.openmc_get_keff.argtypes = [POINTER(c_double*2)]
    _dll.openmc_get_keff.restype = c_int
    _dll.openmc_get_keff.errcheck = _error_handler
    _dll.openmc_plot_geometry.restype = None
    _dll.openmc_run.restype = None
    _dll.openmc_reset.restype = None

    from .nuclide import *
    from .material import *
    from .cell import *
    from .tally import *
