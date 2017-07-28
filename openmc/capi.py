"""Provides bindings to C functions defined by OpenMC shared library.

When the :mod:`openmc` package is imported, the OpenMC shared library is
automatically loaded. Calls to the OpenMC library can then be made, for example:

.. code-block:: python

    openmc.capi.init()
    openmc.capi.run()

"""

from contextlib import contextmanager
from ctypes import CDLL, c_int, c_int32, c_double, c_char_p, POINTER
import sys
from warnings import warn

import numpy as np
from numpy.ctypeslib import as_array
import pkg_resources

_int3 = c_int*3
_double3 = c_double*3
_int_array = POINTER(POINTER(c_int))
_double_array = POINTER(POINTER(c_double))


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


def cell_set_temperature(cell_id, T, instance=None):
    """Set the temperature of a cell

    Parameters
    ----------
    cell_id : int
        ID of the cell
    T : float
        Temperature in K
    instance : int or None
        Which instance of the cell

    """
    _dll.openmc_cell_set_temperature(cell_id, T, instance)



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
    _dll.openmc_find(_double3(*xyz), r_int, uid, instance)
    return (uid.value if uid != 0 else None), instance.value


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


def load_nuclide(name):
    """Load cross section data for a nuclide.

    Parameters
    ----------
    name : str
        Name of nuclide, e.g. 'U235'

    """
    _dll.openmc_load_nuclide(name.encode())


def material_add_nuclide(mat_id, name, density):
    """Add a nuclide to a material.

    Parameters
    ----------
    mat_id : int
        ID of the material
    name : str
        Name of nuclide, e.g. 'U235'
    density : float
        Density in atom/b-cm

    """
    _dll.openmc_material_add_nuclide(mat_id, name.encode(), density)


def material_get_densities(mat_id):
    """Get atom densities in a material.

    Parameters
    ----------
    mat_id : int
        ID of the material

    Returns
    -------
    list of string
        List of nuclide names
    numpy.ndarray
        Array of densities in atom/b-cm

    """
    # Allocate memory for arguments that are written to
    nuclides = POINTER(c_int)()
    densities = POINTER(c_double)()
    n = c_int()

    # Get nuclide names and densities
    _dll.openmc_material_get_densities(mat_id, nuclides, densities, n)

    # Convert to appropriate types and return
    nuclide_list = [nuclide_name(nuclides[i]) for i in range(n.value)]
    density_array = as_array(densities, (n.value,))
    return nuclide_list, density_array


def material_set_density(mat_id, density):
    """Set density of a material.

    Parameters
    ----------
    mat_id : int
        ID of the material
    density : float
        Density in atom/b-cm

    """
    _dll.openmc_material_set_density(mat_id, density)


def material_set_densities(mat_id, nuclides, densities):
    """Set the densities of a list of nuclides in a material

    Parameters
    ----------
    mat_id : int
        ID of the material
    nuclides : iterable of str
        Nuclide names
    densities : iterable of float
        Corresponding densities in atom/b-cm

    """
    # Convert strings to an array of char*
    nucs = (c_char_p * len(nuclides))()
    nucs[:] = [x.encode() for x in nuclides]

    # Get numpy array as a double*
    d = np.asarray(densities)
    dp = d.ctypes.data_as(POINTER(c_double))

    _dll.openmc_material_set_densities(mat_id, len(nuclides), nucs, dp)


def nuclide_name(index):
    """Name of nuclide with given index

    Parameter
    ---------
    index : int
        Index in internal nuclides array

    Returns
    -------
    str
        Name of nuclide

    """
    name = c_char_p()
    _dll.openmc_nuclide_name(index, name)

    # Find blank in name
    i = 0
    while name.value[i:i+1] != b' ':
        i += 1
    return name.value[:i].decode()


def plot_geometry():
    """Plot geometry"""
    _dll.openmc_plot_geometry()


def reset():
    """Reset tallies"""
    _dll.openmc_reset()


def run():
    """Run simulation"""
    _dll.openmc_run()


def tally_results(tally_id):
    """Get tally results array

    Parameters
    ----------
    tally_id : int
        ID of tally

    Returns
    -------
    numpy.ndarray
        Array that exposes the internal tally results array

    """
    data = POINTER(c_double)()
    shape = _int3()
    _dll.openmc_tally_results(tally_id, data, shape)
    return as_array(data, tuple(shape[::-1]))


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

if _available:
    # Set argument/return types
    _dll.openmc_calculate_volumes.restype = None
    _dll.openmc_cell_set_temperature.argtypes = [
        c_int32, c_double, POINTER(c_int32)]
    _dll.openmc_cell_set_temperature.restype = c_int
    _dll.openmc_cell_set_temperature.errcheck = _error_handler
    _dll.openmc_finalize.restype = None
    _dll.openmc_find.argtypes = [
        POINTER(_double3), c_int, POINTER(c_int32), POINTER(c_int32)]
    _dll.openmc_find.restype = c_int
    _dll.openmc_find.errcheck = _error_handler
    _dll.openmc_init.argtypes = [POINTER(c_int)]
    _dll.openmc_init.restype = None
    _dll.openmc_get_keff.argtypes = [POINTER(c_double*2)]
    _dll.openmc_get_keff.restype = c_int
    _dll.openmc_get_keff.errcheck = _error_handler
    _dll.openmc_load_nuclide.argtypes = [c_char_p]
    _dll.openmc_load_nuclide.restype = c_int
    _dll.openmc_load_nuclide.errcheck = _error_handler

    # Material interface
    _dll.openmc_material_add_nuclide.argtypes = [
        c_int32, c_char_p, c_double]
    _dll.openmc_material_add_nuclide.restype = c_int
    _dll.openmc_material_add_nuclide.errcheck = _error_handler
    _dll.openmc_material_get_densities.argtypes = [
        c_int32, _int_array, _double_array, POINTER(c_int)]
    _dll.openmc_material_get_densities.restype = c_int
    _dll.openmc_material_get_densities.errcheck = _error_handler
    _dll.openmc_material_set_density.argtypes = [c_int32, c_double]
    _dll.openmc_material_set_density.restype = c_int
    _dll.openmc_material_set_density.errcheck = _error_handler
    _dll.openmc_material_set_densities.argtypes = [
        c_int32, c_int, POINTER(c_char_p), POINTER(c_double)]
    _dll.openmc_material_set_densities.restype = c_int
    _dll.openmc_material_set_densities.errcheck = _error_handler

    _dll.openmc_nuclide_name.argtypes = [c_int, POINTER(c_char_p)]
    _dll.openmc_nuclide_name.restype = c_int
    _dll.openmc_nuclide_name.errcheck = _error_handler
    _dll.openmc_plot_geometry.restype = None
    _dll.openmc_run.restype = None
    _dll.openmc_reset.restype = None
    _dll.openmc_tally_results.argtypes = [
        c_int32, _double_array, POINTER(_int3)]
    _dll.openmc_tally_results.restype = c_int
    _dll.openmc_tally_results.errcheck = _error_handler
