from ctypes import CDLL, c_int, POINTER, byref, c_double
import sys
from warnings import warn

import numpy as np
from numpy.ctypeslib import ndpointer, as_array

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
        self._dll.openmc_tally_results.argtypes = [
            c_int, POINTER(POINTER(c_double)), ndpointer(
                np.intc, shape=(3,))]
        self._dll.openmc_tally_results.restype = None
        self._dll.openmc_cell_set_temperature.argtypes = [
            c_int, c_double]
        self._dll.openmc_cell_set_temperature.restype = c_int

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

    def tally_results(self, tally_id):
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
        r_p = POINTER(c_double)()
        r_shape = np.zeros(3, np.intc)
        self._dll.openmc_tally_results(tally_id, byref(r_p), r_shape)
        if r_p:
            return as_array(r_p, tuple(r_shape[::-1]))
        else:
            return None

    def cell_set_temperature(self, cell_id, temperature):
        """Set the temperature of a cell"""
        return self._dll.openmc_cell_set_temperature(cell_id, temperature)

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
