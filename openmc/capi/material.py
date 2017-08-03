from collections import Mapping
from ctypes import c_int, c_int32, c_double, c_char_p, POINTER
from weakref import WeakValueDictionary

import numpy as np
from numpy.ctypeslib import as_array

from . import _dll, NuclideView
from .error import _error_handler


__all__ = ['MaterialView', 'materials']

# Material functions
_dll.openmc_get_material.argtypes = [c_int32, POINTER(c_int32)]
_dll.openmc_get_material.restype = c_int
_dll.openmc_get_material.errcheck = _error_handler
_dll.openmc_material_add_nuclide.argtypes = [
    c_int32, c_char_p, c_double]
_dll.openmc_material_add_nuclide.restype = c_int
_dll.openmc_material_add_nuclide.errcheck = _error_handler
_dll.openmc_material_get_id.argtypes = [c_int32, POINTER(c_int32)]
_dll.openmc_material_get_id.restype = c_int
_dll.openmc_material_get_id.errcheck = _error_handler
_dll.openmc_material_get_densities.argtypes = [
    c_int32, POINTER(POINTER(c_int)), POINTER(POINTER(c_double)),
    POINTER(c_int)]
_dll.openmc_material_get_densities.restype = c_int
_dll.openmc_material_get_densities.errcheck = _error_handler
_dll.openmc_material_set_density.argtypes = [c_int32, c_double]
_dll.openmc_material_set_density.restype = c_int
_dll.openmc_material_set_density.errcheck = _error_handler
_dll.openmc_material_set_densities.argtypes = [
    c_int32, c_int, POINTER(c_char_p), POINTER(c_double)]
_dll.openmc_material_set_densities.restype = c_int
_dll.openmc_material_set_densities.errcheck = _error_handler


class MaterialView(object):
    """View of a material.

    This class exposes a material that is stored internally in the OpenMC
    solver. To obtain a view of a material with a given ID, use the
    :data:`openmc.capi.materials` mapping.

    Parameters
    ----------
    index : int
         Index in the `materials` array.

    Attributes
    ----------
    id : int
        ID of the material
    nuclides : list of str
        List of nuclides in the material
    densities : numpy.ndarray
        Array of densities in atom/b-cm

    """
    __instances = WeakValueDictionary()

    def __new__(cls, *args):
        if args not in cls.__instances:
            instance = super().__new__(cls)
            cls.__instances[args] = instance
        return cls.__instances[args]

    def __init__(self, index):
        self._index = index

    @property
    def id(self):
        mat_id = c_int32()
        _dll.openmc_material_get_id(self._index, mat_id)
        return mat_id.value

    @property
    def nuclides(self):
        return self._get_densities()[0]
        return nuclides

    @property
    def densities(self):
        return self._get_densities()[1]

    def _get_densities(self):
        """Get atom densities in a material.

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
        _dll.openmc_material_get_densities(self._index, nuclides, densities, n)

        # Convert to appropriate types and return
        nuclide_list = [NuclideView(nuclides[i]).name for i in range(n.value)]
        density_array = as_array(densities, (n.value,))
        return nuclide_list, density_array

    def add_nuclide(name, density):
        """Add a nuclide to a material.

        Parameters
        ----------
        name : str
            Name of nuclide, e.g. 'U235'
        density : float
            Density in atom/b-cm

        """
        _dll.openmc_material_add_nuclide(self._index, name.encode(), density)

    def set_density(self, density):
        """Set density of a material.

        Parameters
        ----------
        density : float
            Density in atom/b-cm

        """
        _dll.openmc_material_set_density(self._index, density)

    def set_densities(self, nuclides, densities):
        """Set the densities of a list of nuclides in a material

        Parameters
        ----------
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

        _dll.openmc_material_set_densities(self._index, len(nuclides), nucs, dp)


class _MaterialMapping(Mapping):
    def __getitem__(self, key):
        index = c_int32()
        _dll.openmc_get_material(key, index)
        return MaterialView(index.value)

    def __iter__(self):
        for i in range(len(self)):
            yield MaterialView(i + 1).id

    def __len__(self):
        return c_int32.in_dll(_dll, 'n_materials').value

materials = _MaterialMapping()
