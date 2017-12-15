from contextlib import contextmanager
from ctypes import (CDLL, c_int, c_int32, c_int64, c_double, c_char_p,
                    POINTER, Structure)
from warnings import warn

import numpy as np
from numpy.ctypeslib import as_array

from . import _dll
from .error import _error_handler, AllocationError
import openmc.capi


class _Bank(Structure):
    _fields_ = [('wgt', c_double),
                ('xyz', c_double*3),
                ('uvw', c_double*3),
                ('E', c_double),
                ('delayed_group', c_int)]


_dll.openmc_calculate_volumes.restype = None
_dll.openmc_finalize.restype = None
_dll.openmc_find.argtypes = [POINTER(c_double*3), c_int, POINTER(c_int32),
                             POINTER(c_int32)]
_dll.openmc_find.restype = c_int
_dll.openmc_find.errcheck = _error_handler
_dll.openmc_hard_reset.restype = None
_dll.openmc_init.argtypes = [POINTER(c_int)]
_dll.openmc_init.restype = None
_dll.openmc_get_keff.argtypes = [POINTER(c_double*2)]
_dll.openmc_get_keff.restype = c_int
_dll.openmc_get_keff.errcheck = _error_handler
_dll.openmc_next_batch.restype = c_int
_dll.openmc_plot_geometry.restype = None
_dll.openmc_run.restype = None
_dll.openmc_reset.restype = None
_dll.openmc_source_bank.argtypes = [POINTER(POINTER(_Bank)), POINTER(c_int64)]
_dll.openmc_source_bank.restype = c_int
_dll.openmc_source_bank.errcheck = _error_handler
_dll.openmc_simulation_init.restype = None
_dll.openmc_simulation_finalize.restype = None
_dll.openmc_statepoint_write.argtypes = [POINTER(c_char_p)]
_dll.openmc_statepoint_write.restype = None


def calculate_volumes():
    """Run stochastic volume calculation"""
    _dll.openmc_calculate_volumes()


def finalize():
    """Finalize simulation and free memory"""
    _dll.openmc_finalize()


def find_cell(xyz):
    """Find the cell at a given point

    Parameters
    ----------
    xyz : iterable of float
        Cartesian coordinates of position

    Returns
    -------
    openmc.capi.Cell
        Cell containing the point
    int
        If the cell at the given point is repeated in the geometry, this
        indicates which instance it is, i.e., 0 would be the first instance.

    """
    uid = c_int32()
    instance = c_int32()
    _dll.openmc_find((c_double*3)(*xyz), 1, uid, instance)
    return openmc.capi.cells[uid.value], instance.value


def find_material(xyz):
    """Find the material at a given point

    Parameters
    ----------
    xyz : iterable of float
        Cartesian coordinates of position

    Returns
    -------
    openmc.capi.Material or None
        Material containing the point, or None is no material is found

    """
    uid = c_int32()
    instance = c_int32()
    _dll.openmc_find((c_double*3)(*xyz), 2, uid, instance)
    return openmc.capi.materials[uid.value] if uid != 0 else None


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


def iter_batches():
    """Iterator over batches.

    This function returns a generator-iterator that allows Python code to be run
    between batches in an OpenMC simulation. It should be used in conjunction
    with :func:`openmc.capi.simulation_init` and
    :func:`openmc.capi.simulation_finalize`. For example:

    .. code-block:: Python

        with openmc.capi.run_in_memory():
            openmc.capi.simulation_init()
            for _ in openmc.capi.iter_batches():
                # Look at convergence of tallies, for example
                ...
            openmc.capi.simulation_finalize()

    See Also
    --------
    openmc.capi.next_batch

    """
    while True:
        # Run next batch
        retval = next_batch()

        # Provide opportunity for user to perform action between batches
        yield

        # End the iteration
        if retval < 0:
            break


def keff():
    """Return the calculated k-eigenvalue and its standard deviation.

    Returns
    -------
    tuple
        Mean k-eigenvalue and standard deviation of the mean

    """
    n = openmc.capi.num_realizations()
    if n > 3:
        # Use the combined estimator if there are enough realizations
        k = (c_double*2)()
        _dll.openmc_get_keff(k)
        return tuple(k)
    else:
        # Otherwise, return the tracklength estimator
        mean = c_double.in_dll(_dll, 'keff').value
        std_dev = c_double.in_dll(_dll, 'keff_std').value if n > 1 else np.inf
        return (mean, std_dev)


def next_batch():
    """Run next batch."""
    retval = _dll.openmc_next_batch()
    if retval == -3:
        raise AllocationError('Simulation has not been initialized. You must call '
                              'openmc.capi.simulation_init() first.')
    return retval


def plot_geometry():
    """Plot geometry"""
    _dll.openmc_plot_geometry()


def reset():
    """Reset tallies and timers."""
    _dll.openmc_reset()


def run():
    """Run simulation"""
    _dll.openmc_run()


def simulation_init():
    """Initialize simulation"""
    _dll.openmc_simulation_init()


def simulation_finalize():
    """Finalize simulation"""
    _dll.openmc_simulation_finalize()


def source_bank():
    """Return source bank as NumPy array

    Returns
    -------
    numpy.ndarray
        Source sites

    """
    # Get pointer to source bank
    ptr = POINTER(_Bank)()
    n = c_int64()
    _dll.openmc_source_bank(ptr, n)

    # Convert to numpy array with appropriate datatype
    bank_dtype = np.dtype(_Bank)
    return as_array(ptr, (n.value,)).view(bank_dtype)


def statepoint_write(filename=None):
    """Write a statepoint file.

    Parameters
    ----------
    filename : str or None
        Path to the statepoint to write. If None is passed, a default name that
        contains the current batch will be written.

    """
    if filename is not None:
        filename = c_char_p(filename.encode())
    _dll.openmc_statepoint_write(filename)


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
    try:
        yield
    finally:
        finalize()


class _DLLGlobal(object):
    """Data descriptor that exposes global variables from libopenmc."""
    def __init__(self, ctype, name):
        self.ctype = ctype
        self.name = name

    def __get__(self, instance, owner):
        return self.ctype.in_dll(_dll, self.name).value

    def __set__(self, instance, value):
        self.ctype.in_dll(_dll, self.name).value = value


class _FortranObject(object):
    def __repr__(self):
        return "{}[{}]".format(type(self).__name__, self._index)


class _FortranObjectWithID(_FortranObject):
    def __init__(self, uid=None, new=True, index=None):
        # Creating the object has already been handled by __new__. In the
        # initializer, all we do is make sure that the object returned has an ID
        # assigned. If the array index of the object is out of bounds, an
        # OutOfBoundsError will be raised here by virtue of referencing self.id
        self.id
