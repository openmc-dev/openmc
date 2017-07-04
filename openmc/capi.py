from ctypes import CDLL, c_int, POINTER, byref, c_double, c_char_p
import sys
from warnings import warn

import numpy as np
from numpy.ctypeslib import as_array

import pkg_resources

_int3 = c_int*3
_double3 = c_double*3
_double_array = POINTER(POINTER(c_double))


class OpenMCLibrary(object):
    """Provides bindings to C functions defined by OpenMC shared library.

    This class is normally not directly instantiated. Instead, when the
    :mod:`openmc` package is imported, an instance is automatically created with
    the name :data:`openmc.lib`. Calls to the OpenMC can then be made using that
    instance, for example:

    .. code-block:: python

        openmc.lib.init()
        openmc.lib.run()

    """
    def __init__(self, filename):
        self._dll = CDLL(filename)

        # Set argument/return types
        self._dll.openmc_calculate_volumes.restype = None
        self._dll.openmc_cell_set_temperature.argtypes = [
            c_int, c_double]
        self._dll.openmc_cell_set_temperature.restype = c_int
        self._dll.openmc_finalize.restype = None
        self._dll.openmc_find.argtypes = [
            POINTER(_double3), c_int, POINTER(c_int), POINTER(c_int)]
        self._dll.openmc_find.restype = None
        self._dll.openmc_init.argtypes = [POINTER(c_int)]
        self._dll.openmc_init.restype = None
        self._dll.openmc_load_nuclide.argtypes = [c_char_p]
        self._dll.openmc_load_nuclide.restype = c_int
        self._dll.openmc_material_add_nuclide.argtypes = [
            c_int, c_char_p, c_double]
        self._dll.openmc_material_add_nuclide.restype = c_int
        self._dll.openmc_material_get_densities.argtypes = [
            c_int, _double_array]
        self._dll.openmc_material_get_densities.restype = c_int
        self._dll.openmc_material_set_density.argtypes = [c_int, c_double]
        self._dll.openmc_material_set_density.restype = c_int
        self._dll.openmc_plot_geometry.restype = None
        self._dll.openmc_run.restype = None
        self._dll.openmc_reset.restype = None
        self._dll.openmc_set_density.argtypes = [POINTER(_double3), c_double]
        self._dll.openmc_set_density.restype = c_int
        self._dll.openmc_set_temperature.argtypes = [
            POINTER(_double3), c_double]
        self._dll.openmc_set_temperature.restype = c_int
        self._dll.openmc_tally_results.argtypes = [
            c_int, _double_array, POINTER(_int3)]
        self._dll.openmc_tally_results.restype = None

    def calculate_volumes(self):
        """Run stochastic volume calculation"""
        return self._dll.openmc_calculate_volumes()

    def cell_set_temperature(self, cell_id, T):
        """Set the temperature of a cell

        Parameters
        ----------
        cell_id : int
            ID of the cell
        T : float
            Temperature in K

        """
        return self._dll.openmc_cell_set_temperature(cell_id, T)

    def finalize(self):
        """Finalize simulation and free memory"""
        return self._dll.openmc_finalize()

    def find(self, xyz, rtype='cell'):
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
        uid = c_int()
        instance = c_int()
        self._dll.openmc_find(_double3(*xyz), r_int, byref(uid), byref(instance))
        return (uid.value if uid != 0 else None), instance.value

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

    def load_nuclide(self, name):
        """Load cross section data for a nuclide.

        Parameters
        ----------
        name : str
            Name of nuclide, e.g. 'U235'

        Returns
        -------
        int
            Return status (negative if an error occurs).

        """
        return self._dll.openmc_load_nuclide(name.encode())

    def material_add_nuclide(self, mat_id, name, density):
        """Add a nuclide to a material.

        Parameters
        ----------
        mat_id : int
            ID of the material
        name : str
            Name of nuclide, e.g. 'U235'
        density : float
            Density in atom/b-cm

        Returns
        -------
        int
            Return status (negative if an error occurs).

        """
        return self._dll.openmc_material_add_nuclide(
            mat_id, name.encode(), density)

    def material_get_densities(self, mat_id):
        """Get atom densities in a material.

        Parameters
        ----------
        mat_id : int
            ID of the material

        Returns
        -------
        numpy.ndarray
            Array of densities in atom/b-cm

        """
        data = POINTER(c_double)()
        n = self._dll.openmc_material_get_densities(mat_id, byref(data))
        if data:
            return as_array(data, (n,))
        else:
            return None

    def material_set_density(self, mat_id, density):
        """Set density of a material.

        Parameters
        ----------
        mat_id : int
            ID of the material
        density : float
            Density in atom/b-cm

        Returns
        -------
        int
            Return status (negative if an error occurs).

        """
        return self._dll.openmc_material_set_density(mat_id, density)

    def plot_geometry(self):
        """Plot geometry"""
        return self._dll.openmc_plot_geometry()

    def reset(self):
        """Reset tallies"""
        return self._dll.openmc_reset()

    def run(self):
        """Run simulation"""
        return self._dll.openmc_run()

    def set_density(self, xyz, density):
        """Set density at a given point.

        Parameters
        ----------
        xyz : iterable of float
            Cartesian coordinates at position
        density : float
            Density in atom/b-cm

        Returns
        -------
        int
            Return status (negative if an error occurs)

        """
        return self._dll.openmc_set_density(_double3(*xyz), density)

    def set_temperature(self, xyz, T):
        """Set temperature at a given point.

        Parameters
        ----------
        xyz : iterable of float
            Cartesian coordinates at position
        T : float
            Temperature in K

        Returns
        -------
        int
            Return status (negative if an error occurs)

        """
        return self._dll.openmc_set_temperature(_double3(*xyz), T)

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
        data = POINTER(c_double)()
        shape = _int3()
        self._dll.openmc_tally_results(tally_id, byref(data), byref(shape))
        if data:
            return as_array(data, tuple(shape[::-1]))
        else:
            return None

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
    lib = OpenMCLibrary(filename)
except OSError:
    warn("OpenMC shared library is not available from the Python API. This "
         "means you will not be able to use openmc.lib to make in-memory "
         "calls to OpenMC.")
    lib = None
