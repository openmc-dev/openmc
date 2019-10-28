from collections.abc import Mapping
from ctypes import c_int, c_int32, c_double, c_char_p, POINTER, c_size_t
from weakref import WeakValueDictionary

import numpy as np
from numpy.ctypeslib import as_array

from openmc.exceptions import AllocationError, InvalidIDError, OpenMCError
from . import _dll, Nuclide
from .core import _FortranObjectWithID
from .error import _error_handler


__all__ = ['Material', 'materials']

# Material functions
_dll.openmc_extend_materials.argtypes = [c_int32, POINTER(c_int32), POINTER(c_int32)]
_dll.openmc_extend_materials.restype = c_int
_dll.openmc_extend_materials.errcheck = _error_handler
_dll.openmc_get_material_index.argtypes = [c_int32, POINTER(c_int32)]
_dll.openmc_get_material_index.restype = c_int
_dll.openmc_get_material_index.errcheck = _error_handler
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
_dll.openmc_material_get_density.argtypes = [c_int32, POINTER(c_double)]
_dll.openmc_material_get_density.restype = c_int
_dll.openmc_material_get_density.errcheck = _error_handler
_dll.openmc_material_get_volume.argtypes = [c_int32, POINTER(c_double)]
_dll.openmc_material_get_volume.restype = c_int
_dll.openmc_material_get_volume.errcheck = _error_handler
_dll.openmc_material_set_density.argtypes = [c_int32, c_double, c_char_p]
_dll.openmc_material_set_density.restype = c_int
_dll.openmc_material_set_density.errcheck = _error_handler
_dll.openmc_material_set_densities.argtypes = [
    c_int32, c_int, POINTER(c_char_p), POINTER(c_double)]
_dll.openmc_material_set_densities.restype = c_int
_dll.openmc_material_set_densities.errcheck = _error_handler
_dll.openmc_material_set_id.argtypes = [c_int32, c_int32]
_dll.openmc_material_set_id.restype = c_int
_dll.openmc_material_set_id.errcheck = _error_handler
_dll.openmc_material_get_name.argtypes = [c_int32, POINTER(c_char_p)]
_dll.openmc_material_get_name.restype = c_int
_dll.openmc_material_get_name.errcheck = _error_handler
_dll.openmc_material_set_name.argtypes = [c_int32, c_char_p]
_dll.openmc_material_set_name.restype = c_int
_dll.openmc_material_set_name.errcheck = _error_handler
_dll.openmc_material_set_volume.argtypes = [c_int32, c_double]
_dll.openmc_material_set_volume.restype = c_int
_dll.openmc_material_set_volume.errcheck = _error_handler
_dll.n_materials.argtypes = []
_dll.n_materials.restype = c_size_t


class Material(_FortranObjectWithID):
    """Material stored internally.

    This class exposes a material that is stored internally in the OpenMC
    library. To obtain a view of a material with a given ID, use the
    :data:`openmc.lib.materials` mapping.

    Parameters
    ----------
    uid : int or None
        Unique ID of the tally
    new : bool
        When `index` is None, this argument controls whether a new object is
        created or a view to an existing object is returned.
    index : int or None
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

    def __new__(cls, uid=None, new=True, index=None):
        mapping = materials
        if index is None:
            if new:
                # Determine ID to assign
                if uid is None:
                    uid = max(mapping, default=0) + 1
                else:
                    if uid in mapping:
                        raise AllocationError('A material with ID={} has already '
                                              'been allocated.'.format(uid))

                index = c_int32()
                _dll.openmc_extend_materials(1, index, None)
                index = index.value
            else:
                index = mapping[uid]._index
        elif index == -1:
            # Special value indicates void material
            return None

        if index not in cls.__instances:
            instance = super(Material, cls).__new__(cls)
            instance._index = index
            if uid is not None:
                instance.id = uid
            cls.__instances[index] = instance

        return cls.__instances[index]

    @property
    def id(self):
        mat_id = c_int32()
        _dll.openmc_material_get_id(self._index, mat_id)
        return mat_id.value

    @id.setter
    def id(self, mat_id):
        _dll.openmc_material_set_id(self._index, mat_id)

    @property
    def name(self):
        name = c_char_p()
        _dll.openmc_material_get_name(self._index, name)
        return name.value.decode()

    @name.setter
    def name(self, name):
        name_ptr = c_char_p(name.encode())
        _dll.openmc_material_set_name(self._index, name_ptr)

    @property
    def volume(self):
        volume = c_double()
        try:
            _dll.openmc_material_get_volume(self._index, volume)
        except OpenMCError:
            return None
        return volume.value

    @volume.setter
    def volume(self, volume):
        _dll.openmc_material_set_volume(self._index, volume)

    @property
    def nuclides(self):
        return self._get_densities()[0]
        return nuclides

    @property
    def density(self):
      density = c_double()
      try:
          _dll.openmc_material_get_density(self._index, density)
      except OpenMCError:
          return None
      return density.value

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
        nuclide_list = [Nuclide(nuclides[i]).name for i in range(n.value)]
        density_array = as_array(densities, (n.value,))
        return nuclide_list, density_array

    def add_nuclide(self, name, density):
        """Add a nuclide to a material.

        Parameters
        ----------
        name : str
            Name of nuclide, e.g. 'U235'
        density : float
            Density in atom/b-cm

        """
        _dll.openmc_material_add_nuclide(self._index, name.encode(), density)

    def set_density(self, density, units='atom/b-cm'):
        """Set density of a material.

        Parameters
        ----------
        density : float
            Density
        units : {'atom/b-cm', 'g/cm3'}
            Units for density

        """
        _dll.openmc_material_set_density(self._index, density, units.encode())

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
        try:
            _dll.openmc_get_material_index(key, index)
        except (AllocationError, InvalidIDError) as e:
            # __contains__ expects a KeyError to work correctly
            raise KeyError(str(e))
        return Material(index=index.value)

    def __iter__(self):
        for i in range(len(self)):
            yield Material(index=i).id

    def __len__(self):
        return _dll.n_materials()

    def __repr__(self):
        return repr(dict(self))

materials = _MaterialMapping()
