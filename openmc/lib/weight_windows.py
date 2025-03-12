from collections.abc import Mapping
from ctypes import c_double, c_int, c_int32, c_char_p, c_size_t, POINTER
from weakref import WeakValueDictionary

import numpy as np
from numpy.ctypeslib import as_array

from openmc import ParticleType
from openmc.exceptions import AllocationError, InvalidIDError
from . import _dll
from .core import _FortranObjectWithID
from .error import _error_handler
from .filter import EnergyFilter, MeshFilter, ParticleFilter
from .mesh import _get_mesh
from .mesh import meshes


__all__ = ['WeightWindows', 'weight_windows']

_dll.openmc_extend_weight_windows.argtypes = [c_int32, POINTER(c_int32), POINTER(c_int32)]

_dll.openmc_weight_windows_update_magic.argtypes = 2*[c_int32] + [c_char_p] + 2*[c_double]
_dll.openmc_weight_windows_update_magic.restype = c_int
_dll.openmc_weight_windows_update_magic.errcheck = _error_handler

_dll.openmc_weight_windows_size.restype = c_size_t

_dll.openmc_get_weight_windows_index.argtypes = [c_int32, POINTER(c_int32)]
_dll.openmc_get_weight_windows_index.restype = c_int
_dll.openmc_get_weight_windows_index.errcheck = _error_handler

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

_dll.openmc_weight_windows_set_particle.argtypes = [c_int32, c_int]
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

_dll.openmc_weight_windows_get_survival_ratio.argtypes = [c_int32, POINTER(c_double)]
_dll.openmc_weight_windows_get_survival_ratio.restype = c_int
_dll.openmc_weight_windows_get_survival_ratio.errcheck = _error_handler

_dll.openmc_weight_windows_set_survival_ratio.argtypes = [c_int32, c_double]
_dll.openmc_weight_windows_set_survival_ratio.restype = c_int
_dll.openmc_weight_windows_set_survival_ratio.errcheck = _error_handler

_dll.openmc_weight_windows_get_max_lower_bound_ratio.argtypes = [c_int32, POINTER(c_double)]
_dll.openmc_weight_windows_get_max_lower_bound_ratio.restype = c_int
_dll.openmc_weight_windows_get_max_lower_bound_ratio.errcheck = _error_handler

_dll.openmc_weight_windows_set_max_lower_bound_ratio.argtypes = [c_int32, c_double]
_dll.openmc_weight_windows_set_max_lower_bound_ratio.restype = c_int
_dll.openmc_weight_windows_set_max_lower_bound_ratio.errcheck = _error_handler

_dll.openmc_weight_windows_get_weight_cutoff.argtypes = [c_int32, POINTER(c_double)]
_dll.openmc_weight_windows_get_weight_cutoff.restype = c_int
_dll.openmc_weight_windows_get_weight_cutoff.errcheck = _error_handler

_dll.openmc_weight_windows_set_weight_cutoff.argtypes = [c_int32, c_double]
_dll.openmc_weight_windows_set_weight_cutoff.restype = c_int
_dll.openmc_weight_windows_set_weight_cutoff.errcheck = _error_handler

_dll.openmc_weight_windows_get_max_split.argtypes = [c_int32, POINTER(c_int)]
_dll.openmc_weight_windows_get_max_split.restype = c_int
_dll.openmc_weight_windows_get_max_split.errcheck = _error_handler

_dll.openmc_weight_windows_set_max_split.argtypes = [c_int32, c_int]
_dll.openmc_weight_windows_set_max_split.restype = c_int
_dll.openmc_weight_windows_set_max_split.errcheck = _error_handler


class WeightWindows(_FortranObjectWithID):
    """WeightWindows stored internally.

    This class exposes a weight windows object that is stored internally in the
    OpenMC library. To obtain a view of a weight windows object with a given ID,
    use the :data:`openmc.lib.weight_windows` mapping.

    .. versionadded:: 0.14.0

    Parameters
    ----------
    id : int or None
        Unique ID of the weight windows
    new : bool
        When `index` is None, this argument controls whether a new object is
        created or a view of an existing object is returned.
    index : int or None
        Index in the `weight_windows` array.

    Attributes
    ----------
    id : int
        ID of the weight windows object
    mesh : openmc.lib.Mesh
        Mesh used for the weight windows
    particle : openmc.ParticleType
        The particle type to which these weight windows apply
    energy_bounds : numpy.ndarray
        The energy bounds for the weight windows
    bounds : numpy.ndarray
        The weight window bounds
    """
    __instances = WeakValueDictionary()

    def __new__(cls, id=None, new=True, index=None):
        mapping = weight_windows

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
        _dll.openmc_weight_windows_set_mesh(
            weight_windows[self.id]._index, meshes[mesh.id]._index)

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
        _dll.openmc_weight_windows_set_energy_bounds(
            self._index, e_bounds_ptr, e_bounds_arr.size)

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
        _dll.openmc_weight_windows_set_particle(self._index, int(p))

    @property
    def bounds(self):
        upper = POINTER(c_double)()
        lower = POINTER(c_double)()
        size = c_size_t()
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

        _dll.openmc_weight_windows_set_bounds(self._index, lower_p, upper_p, lower.size)

    @property
    def survival_ratio(self):
        ratio = c_double()
        _dll.openmc_weight_windows_get_survival_ratio(self._index, ratio)
        return ratio.value

    @survival_ratio.setter
    def survival_ratio(self, ratio):
        _dll.openmc_weight_windows_set_survival_ratio(self._index, ratio)

    @property
    def max_lower_bound_ratio(self):
        lb_ratio = c_double()
        _dll.openmc_weight_windows_get_max_lower_bound_ratio(self._index, lb_ratio)
        return lb_ratio.value

    @max_lower_bound_ratio.setter
    def max_lower_bound_ratio(self, lb_ratio):
        _dll.openmc_weight_windows_set_max_lower_bound_ratio(self._index, lb_ratio)

    @property
    def weight_cutoff(self):
        cutoff = c_double()
        _dll.openmc_weight_windows_get_weight_cutoff(self._index, cutoff)
        return cutoff.value

    @weight_cutoff.setter
    def weight_cutoff(self, cutoff):
        _dll.openmc_weight_windows_set_weight_cutoff(self._index, cutoff)

    @property
    def max_split(self):
        max_split = c_int()
        _dll.openmc_weight_windows_get_max_split(self._index, max_split)
        return max_split.value

    @max_split.setter
    def max_split(self, max_split):
        _dll.openmc_weight_windows_set_max_split(self._index, max_split)

    def update_magic(self, tally, value='mean', threshold=1.0, ratio=5.0):
        """Update weight window values using the MAGIC method

        Reference: https://inis.iaea.org/records/231pm-zzy35

        Parameters
        ----------
        tally : openmc.lib.Tally object
            The tally used to update weight window information
        value : str
            Value type used to generate weight windows. One of {'mean', 'rel_err'}.
        threshold : float
            Threshold for relative error of results used to generate weight window bounds
        ratio : float
            Ratio of the lower to upper weight window bounds

        """
        _dll.openmc_weight_windows_update_magic(self._index,
                                                tally._index,
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
        particle : openmc.ParticleType or str, optional
            The particle type to use for the WeightWindows instance. Should be
            specified as an instance of ParticleType or as a string with a value of
            'neutron' or 'photon'.

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
        ValueError
            If the specified particle is not included in the bins of the ParticleFilter
            of the tally.
        ValueError
            If the tally does not have a MeshFilter.
        """
        # do some checks on particle value
        if not isinstance(particle, (ParticleType, str)):
            raise ValueError(f"Parameter 'particle' must be {ParticleType} or one of ('neutron', 'photon').")

        # convert particle type if needed
        if isinstance(particle, str):
            particle = ParticleType.from_string(particle)

        if particle not in (ParticleType.NEUTRON, ParticleType.PHOTON):
            raise ValueError('Weight windows can only be applied for neutrons or photons')

        try:
            particle_filter = tally.find_filter(ParticleFilter)
        except ValueError:
            particle_filter = None

        # ensure that the tally won't filter out the specified particle
        if particle_filter is not None and particle not in particle_filter.bins:
            raise ValueError(f'Specified tally for weight windows (Tally {tally.id})'
                             f' does not track the requested particle: "{particle}"')

        # tally must have a mesh filter
        mesh_filter = tally.find_filter(MeshFilter)

        # create a new weight windows instance
        out = cls()

        # set mesh and particle
        out.mesh = mesh_filter.mesh
        out.particle = particle

        # set energy bounds if needed
        try:
            energy_filter = tally.find_filter(EnergyFilter)
        except ValueError:
            energy_filter = None

        if energy_filter is not None:
            out.energy_bounds = energy_filter.bins

        return out


class _WeightWindowsMapping(Mapping):
    def __getitem__(self, key):
        index = c_int32()
        try:
            _dll.openmc_get_weight_windows_index(key, index)
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

weight_windows = _WeightWindowsMapping()
