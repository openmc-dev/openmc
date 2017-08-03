from contextlib import contextmanager
from ctypes import CDLL, c_int, c_int32, c_double, POINTER
from warnings import warn

from . import _dll
from .error import _error_handler


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
_dll.openmc_plot_geometry.restype = None
_dll.openmc_run.restype = None
_dll.openmc_reset.restype = None


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
    int
        ID of the cell.
    int
        If the cell at the given point is repeated in the geometry, this
        indicates which instance it is, i.e., 0 would be the first instance.

    """
    uid = c_int32()
    instance = c_int32()
    _dll.openmc_find((c_double*3)(*xyz), 1, uid, instance)
    return uid.value, instance.value


def find_material(xyz):
    """Find the material at a given point

    Parameters
    ----------
    xyz : iterable of float
        Cartesian coordinates of position

    Returns
    -------
    int or None
        ID of the material or None is no material is found

    """
    uid = c_int32()
    instance = c_int32()
    _dll.openmc_find((c_double*3)(*xyz), 2, uid, instance)
    return uid.value if uid != 0 else None


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
