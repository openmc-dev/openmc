from collections.abc import Mapping
from ctypes import c_double, c_int, c_int32, c_char_p, c_size_t, POINTER
from weakref import WeakValueDictionary

import numpy as np
from numpy.ctypeslib import as_array

from openmc.exceptions import AllocationError, InvalidIDError
from . import _dll
from .core import _FortranObjectWithID
from .error import _error_handler
from .filter import EnergyFilter, MeshFilter, ParticleFilter
from .mesh import _get_mesh
from .mesh import meshes
from .particle import ParticleType


__all__ = ['WeightWindows', 'weight_windows_map']

_dll.openmc_extend_weight_windows.argtypes = [c_int32, POINTER(c_int32), POINTER(c_int32)]

_dll.openmc_update_weight_windows_magic.argtypes = 2*[c_int32] + [c_char_p] + 2*[c_double]
_dll.openmc_update_weight_windows_magic.restype = c_int
_dll.openmc_update_weight_windows_magic.errcheck = _error_handler

_dll.openmc_weight_windows_size.restype = c_size_t

_dll.openmc_weight_windows_get_index.argtypes = [c_int32, POINTER(c_int32)]
_dll.openmc_weight_windows_get_index.restype = c_int
_dll.openmc_weight_windows_get_index.errcheck = _error_handler

_dll.openmc_weight_windows_get_id.argtypes = [c_int32, POINTER(c_int32)]
_dll.openmc_weight_windows_get_id.restype = c_int
_dll.openmc_weight_windows_get_id.errcheck = _error_handler

_dll.openmc_weight_windows_set_id.argtypes = [c_int32, c_int32]
_dll.openmc_weight_windows_set_id.restype = c_int
_dll.openmc_weight_windows_set_id.errcheck = _error_handler

_dll.openmc_weight_windows_set_mesh.argtypes = [c_int32, c_int32]
_dll.openmc_weight_windows_set_mesh.restype = c_int
_dll.openmc_weight_windows_set_mesh.errcheck = _error_handler

_dll.openmc_weight_windows_get_mesh.argtypes = [c_int32, POINTER(c_int32)]
_dll.openmc_weight_windows_get_mesh.restype = c_int
_dll.openmc_weight_windows_get_mesh.errcheck = _error_handler

_dll.openmc_weight_windows_set_energy_bounds.argtypes = [c_int32, POINTER(c_double), c_size_t]
_dll.openmc_weight_windows_set_energy_bounds.restype = c_int
_dll.openmc_weight_windows_set_energy_bounds.errcheck = _error_handler

_dll.openmc_weight_windows_get_energy_bounds.argtypes = [c_int32, POINTER(POINTER(c_double)), POINTER(c_size_t)]
_dll.openmc_weight_windows_get_energy_bounds.restype = c_int
_dll.openmc_weight_windows_get_energy_bounds.errcheck = _error_handler

_dll.openmc_weight_windows_set_particle.argtypes = [c_int32, c_char_p]
_dll.openmc_weight_windows_set_particle.restype = c_int
_dll.openmc_weight_windows_set_particle.errcheck = _error_handler

_dll.openmc_weight_windows_get_particle.argtypes = [c_int32, POINTER(c_int)]
_dll.openmc_weight_windows_get_particle.restype = c_int
_dll.openmc_weight_windows_get_particle.errcheck = _error_handler

_dll.openmc_weight_windows_set_bounds.argtypes = [c_int32, POINTER(c_double), POINTER(c_double), c_size_t]
_dll.openmc_weight_windows_set_bounds.restype = c_int
_dll.openmc_weight_windows_set_bounds.errcheck = _error_handler

_dll.openmc_weight_windows_get_bounds.argtypes = [c_int32, POINTER(POINTER(c_double)), POINTER(POINTER(c_double)), POINTER(c_size_t)]
_dll.openmc_weight_windows_get_bounds.restype = c_int
_dll.openmc_weight_windows_get_bounds.errcheck = _error_handler


class WeightWindows(_FortranObjectWithID):
    """WeightWindows stored internally.

    Parameters
    ----------
    id : int or None
        Unique ID of the weight windows
    new : bool
        When `index` is None, this argument controls whether a new object is
        created or a view of an existing object is returned.
    index : int or None
        Index in the `weight_windows` array.

    """

    __instances = WeakValueDictionary()

    def __new__(cls, id=None, new=True, index=None):
        mapping = weight_windows_map

        if index is None:
            if new:
                # Determine ID to assign
                if id is None:
                    id = max(mapping, default=0) + 1
                else:
                    if id in mapping:
                        raise AllocationError(f'A weight windows object with ID={id} '
                                              'has already been allocated.')

                index = c_int32()
                _dll.openmc_extend_weight_windows(1, index, None)
                index = index.value
            else:
                index = mapping[id]._index

        if index not in cls.__instances:
            instance = super().__new__(cls)
            instance._index = index
            if id is not None:
                instance.id = id
            cls.__instances[index] = instance

        return cls.__instances[index]

    @property
    def id(self):
        ww_id = c_int32()
        _dll.openmc_weight_windows_get_id(self._index, ww_id)
        return ww_id.value

    @id.setter
    def id(self, ww_id):
        _dll.openmc_weight_windows_set_id(self._index, ww_id)

    @property
    def mesh(self):
        mesh_idx = c_int32()
        _dll.openmc_weight_windows_get_mesh(self._index, mesh_idx)
        return _get_mesh(mesh_idx.value)

    @mesh.setter
    def mesh(self, mesh):
        _dll.openmc_weight_windows_set_mesh(weight_windows_map[self.id]._index, meshes[mesh.id]._index)

    @property
    def energy_bounds(self):
        data = POINTER(c_double)()
        n = c_size_t()
        _dll.openmc_weight_windows_get_energy_bounds(self._index, data, n)
        return as_array(data, (n.value,))

    @energy_bounds.setter
    def energy_bounds(self, e_bounds):
        e_bounds_arr = np.asarray(e_bounds, dtype=float)
        e_bounds_ptr = e_bounds_arr.ctypes.data_as(POINTER(c_double))
        _dll.openmc_weight_windows_set_energy_bounds(self._index, e_bounds_ptr, e_bounds_arr.size)

    @property
    def particle(self):
        val = c_int()
        _dll.openmc_weight_windows_get_particle(self._index, val)
        return ParticleType(val.value)

    @particle.setter
    def particle(self, p):
        if isinstance(p, str):
            p = ParticleType.from_string(p)
        else:
            p = ParticleType(p)
        val = c_char_p(str(p).encode())
        _dll.openmc_weight_windows_set_particle(self._index, val)

    @property
    def bounds(self):
        upper = POINTER(c_double)()
        lower = POINTER(c_double)()
        size = c_size_t()
        print(lower)
        _dll.openmc_weight_windows_get_bounds(self._index, lower, upper, size)
        lower_arr = as_array(lower, (size.value,))
        upper_arr = as_array(upper, (size.value,))
        return (lower_arr, upper_arr)

    @bounds.setter
    def bounds(self, bounds):
        lower = np.asarray(bounds[0])
        upper = np.asarray(bounds[1])

        lower_p = lower.ctypes.data_as(POINTER(c_double))
        upper_p = upper.ctypes.data_as(POINTER(c_double))

        _dll.openmc_weight_windows_set_bounds(self._index, lower_p, upper_p, lower_p.size)

    def update_weight_windows_magic(self, tally, value='mean', threshold=1.0, ratio=5.0):
        """Update weight window values using tally information

        Parameters
        ----------
        tally : openmc.lib.Tally object
            The tally used to update weight window information
        value : str
            Value type used to generate weight windows. One of {'mean', 'rel_err', 'std_dev}.
            (default is 'mean')
        threshold : float
            Threshold for relative error of results used to generate weight window bounds
        ratio : float
            Ratio of the lower to upper weight window bounds

        """
        _dll.openmc_update_weight_windows_magic(tally._index,
                                                self._index,
                                                c_char_p(value.encode()),
                                                threshold,
                                                ratio)

    @classmethod
    def from_tally(cls, tally, particle=ParticleType.NEUTRON):
        """Create an instance of the WeightWindows class based on the specified tally.

        Parameters
        ----------
        tally : openmc.lib.Tally
            The tally used to create the WeightWindows instance.
        particle : openmc.lib.ParticleType or str, optional
            The particle type to use for the WeightWindows instance. Should be
            specified as an instance of ParticleType or as a string with a value of
            'neutron' or 'photon'. Defaults to ParticleType.NEUTRON.

        Returns
        -------
        WeightWindows
            The WeightWindows instance created from the specified tally.

        Raises
        ------
        ValueError
            If the particle parameter is not an instance of ParticleType or a string.
        ValueError
            If the particle parameter is not a valid particle type (i.e., not 'neutron'
            or 'photon').
        RuntimeError
            If the specified particle is not included in the bins of the ParticleFilter
            of the tally.
        RuntimeError
            If the tally does not have a MeshFilter.
        """
        # do some checks on particle value
        if not isinstance(particle, (ParticleType, str)):
            raise ValueError(f"Parameter 'particle' must be {ParticleType} or one of ('neutron', 'photon').")

        # convert particle type if needed
        if isinstance(particle, str):
            particle = ParticleType.from_string(particle)

        if particle not in (ParticleType.NEUTRON, ParticleType.PHOTON):
            raise ValueError(f'Weight windows can only be applied for neutrons or photons')

        particle_filter = tally.find_filter(ParticleFilter)
        # ensure that the tally won't filter out the specified particle
        if particle_filter is not None and particle not in particle_filter.bins:
            raise RuntimeError(f'Specified tally for weight windows (Tally {tally.id})'
                               f' does not track the reqeusted particle: "{particle}"')

        # tally has to have a mesh filter
        mesh_filter = tally.find_filter(MeshFilter)
        if mesh_filter is None:
            raise RuntimeError(f'No mesh filter found on tally {tally.id}')

        # create a new weight windows instance
        out = cls()

        # set mesh and particle
        out.mesh = mesh_filter.mesh
        out.particle = particle

        # set energy bounds if needed
        energy_filter = tally.find_filter(EnergyFilter)
        if energy_filter is not None:
            out.energy_bounds = energy_filter.bins

        return out


class _WeightWindowsMapping(Mapping):
    def __getitem__(self, key):
        index = c_int32()
        try:
            _dll.openmc_weight_windows_get_index(key, index)
        except (AllocationError, InvalidIDError) as e:
            raise KeyError(str(e))
        return WeightWindows(index=index.value)

    def __iter__(self):
        for i in range(len(self)):
            yield WeightWindows(index=i).id

    def __len__(self):
        return _dll.openmc_weight_windows_size()

    def __repr__(self):
        return repr(dict(self))

    def __delitem__(self):
        raise NotImplementedError("WeightWindows object remove not implemented")

weight_windows_map = _WeightWindowsMapping()