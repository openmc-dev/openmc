from ctypes import CDLL, c_int, POINTER, byref
import sys
from warnings import warn

import pkg_resources


class _OpenMCLibrary(object):
    def __init__(self, filename):
        self._dll = CDLL(filename)

        # Set argument/return types
        self._dll.openmc_init.argtypes = [POINTER(c_int)]
        self._dll.openmc_init.restype = None
        self._dll.openmc_run.restype = None
        self._dll.openmc_plot_geometry.restype = None
        self._dll.openmc_calculate_volumes.restype = None
        self._dll.openmc_finalize.restype = None
        self._dll.openmc_reset.restype = None

    def init(self, intracomm=None):
        """Initialize OpenMC

        Parameters
        ----------
        intracomm : int or None
            MPI intracommunicator

        """
        if intracomm is not None:
            # If an mpi4py communicator was passed, convert it to an integer to
            # be passed to openmc_init
            try:
                intracomm = intracomm.py2f()
            except AttributeError:
                pass
            return self._dll.openmc_init(byref(c_int(intracomm)))
        else:
            return self._dll.openmc_init(None)

    def run(self):
        """Run simulation"""
        return self._dll.openmc_run()

    def plot_geometry(self):
        """Plot geometry"""
        return self._dll.openmc_plot_geometry()

    def calculate_volumes(self):
        """Run stochastic volume calculation"""
        return self._dll.openmc_calculate_volumes()

    def finalize(self):
        """Finalize simulation and free memory"""
        return self._dll.openmc_finalize()

    def reset(self):
        """Reset tallies"""
        return self._dll.openmc_reset()

    def __getattr__(self, key):
        # Fall-back for other functions that may be available from library
        try:
            return getattr(self._dll, 'openmc_{}'.format(key))
        except AttributeError:
            raise AttributeError("OpenMC library doesn't have a '{}' function"
                                 .format(key))


# Determine shared-library suffix
if sys.platform == 'darwin':
    suffix = 'dylib'
else:
    suffix = 'so'

# Open shared library
filename = pkg_resources.resource_filename(
    __name__, '_libopenmc.{}'.format(suffix))
try:
    lib = _OpenMCLibrary(filename)
except OSError:
    warn("OpenMC shared library is not available from the Python API. This "
         "means you will not be able to use openmc.lib to make in-memory "
         "calls to OpenMC.")
    lib = None
