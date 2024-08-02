from collections.abc import Mapping
from ctypes import c_int, c_int32, c_double, c_char_p, POINTER, \
    create_string_buffer, c_size_t
from weakref import WeakValueDictionary

import numpy as np
from numpy.ctypeslib import as_array

from openmc.exceptions import AllocationError, InvalidIDError
from openmc.data.function import INTERPOLATION_SCHEME
from openmc import ParticleType
from . import _dll
from .core import _FortranObjectWithID
from .error import _error_handler
from .material import Material
from .mesh import _get_mesh


__all__ = [
    'Filter', 'AzimuthalFilter', 'CellFilter', 'CellbornFilter', 'CellfromFilter',
    'CellInstanceFilter', 'CollisionFilter', 'DistribcellFilter', 'DelayedGroupFilter',
    'EnergyFilter', 'EnergyoutFilter', 'EnergyFunctionFilter', 'LegendreFilter',
    'MaterialFilter', 'MaterialFromFilter', 'MeshFilter', 'MeshBornFilter',
    'MeshSurfaceFilter', 'MuFilter', 'ParticleFilter',
    'PolarFilter', 'SphericalHarmonicsFilter', 'SpatialLegendreFilter', 'SurfaceFilter',
    'UniverseFilter', 'ZernikeFilter', 'ZernikeRadialFilter', 'filters'
]

# Tally functions
_dll.openmc_cell_filter_get_bins.argtypes = [
    c_int32, POINTER(POINTER(c_int32)), POINTER(c_int32)]
_dll.openmc_cell_filter_get_bins.restype = c_int
_dll.openmc_cell_filter_get_bins.errcheck = _error_handler
_dll.openmc_energy_filter_get_bins.argtypes = [
    c_int32, POINTER(POINTER(c_double)), POINTER(c_size_t)]
_dll.openmc_energy_filter_get_bins.restype = c_int
_dll.openmc_energy_filter_get_bins.errcheck = _error_handler
_dll.openmc_energy_filter_set_bins.argtypes = [c_int32, c_size_t, POINTER(c_double)]
_dll.openmc_energy_filter_set_bins.restype = c_int
_dll.openmc_energy_filter_set_bins.errcheck = _error_handler
_dll.openmc_energyfunc_filter_set_data.restype = c_int
_dll.openmc_energyfunc_filter_set_data.errcheck = _error_handler
_dll.openmc_energyfunc_filter_set_data.argtypes = [
    c_int32, c_size_t, POINTER(c_double), POINTER(c_double)]
_dll.openmc_energyfunc_filter_get_energy.resttpe = c_int
_dll.openmc_energyfunc_filter_get_energy.errcheck = _error_handler
_dll.openmc_energyfunc_filter_get_energy.argtypes = [
    c_int32, POINTER(c_size_t), POINTER(POINTER(c_double))]
_dll.openmc_energyfunc_filter_get_y.resttpe = c_int
_dll.openmc_energyfunc_filter_get_y.errcheck = _error_handler
_dll.openmc_energyfunc_filter_get_y.argtypes = [
    c_int32, POINTER(c_size_t), POINTER(POINTER(c_double))]
_dll.openmc_energyfunc_filter_get_interpolation.resttpe = c_int
_dll.openmc_energyfunc_filter_get_interpolation.errcheck = _error_handler
_dll.openmc_energyfunc_filter_get_interpolation.argtypes = [c_int32, POINTER(c_int)]
_dll.openmc_energyfunc_filter_set_interpolation.resttpe = c_int
_dll.openmc_energyfunc_filter_set_interpolation.errcheck = _error_handler
_dll.openmc_energyfunc_filter_set_interpolation.argtypes = [c_int32, c_char_p]
_dll.openmc_filter_get_id.argtypes = [c_int32, POINTER(c_int32)]
_dll.openmc_filter_get_id.restype = c_int
_dll.openmc_filter_get_id.errcheck = _error_handler
_dll.openmc_filter_get_num_bins.argtypes = [c_int32, POINTER(c_int)]
_dll.openmc_filter_get_num_bins.restype = c_int
_dll.openmc_filter_get_num_bins.errchck = _error_handler
_dll.openmc_filter_get_type.argtypes = [c_int32, c_char_p]
_dll.openmc_filter_get_type.restype = c_int
_dll.openmc_filter_get_type.errcheck = _error_handler
_dll.openmc_filter_set_id.argtypes = [c_int32, c_int32]
_dll.openmc_filter_set_id.restype = c_int
_dll.openmc_filter_set_id.errcheck = _error_handler
_dll.openmc_get_filter_index.argtypes = [c_int32, POINTER(c_int32)]
_dll.openmc_get_filter_index.restype = c_int
_dll.openmc_get_filter_index.errcheck = _error_handler
_dll.openmc_legendre_filter_get_order.argtypes = [c_int32, POINTER(c_int)]
_dll.openmc_legendre_filter_get_order.restype = c_int
_dll.openmc_legendre_filter_get_order.errcheck = _error_handler
_dll.openmc_legendre_filter_set_order.argtypes = [c_int32, c_int]
_dll.openmc_legendre_filter_set_order.restype = c_int
_dll.openmc_legendre_filter_set_order.errcheck = _error_handler
_dll.openmc_material_filter_get_bins.argtypes = [
    c_int32, POINTER(POINTER(c_int32)), POINTER(c_size_t)]
_dll.openmc_material_filter_get_bins.restype = c_int
_dll.openmc_material_filter_get_bins.errcheck = _error_handler
_dll.openmc_material_filter_set_bins.argtypes = [c_int32, c_size_t, POINTER(c_int32)]
_dll.openmc_material_filter_set_bins.restype = c_int
_dll.openmc_material_filter_set_bins.errcheck = _error_handler
_dll.openmc_mesh_filter_get_mesh.argtypes = [c_int32, POINTER(c_int32)]
_dll.openmc_mesh_filter_get_mesh.restype = c_int
_dll.openmc_mesh_filter_get_mesh.errcheck = _error_handler
_dll.openmc_mesh_filter_set_mesh.argtypes = [c_int32, c_int32]
_dll.openmc_mesh_filter_set_mesh.restype = c_int
_dll.openmc_mesh_filter_set_mesh.errcheck = _error_handler
_dll.openmc_mesh_filter_get_translation.argtypes = [c_int32, POINTER(c_double*3)]
_dll.openmc_mesh_filter_get_translation.restype = c_int
_dll.openmc_mesh_filter_get_translation.errcheck = _error_handler
_dll.openmc_mesh_filter_set_translation.argtypes = [c_int32, POINTER(c_double*3)]
_dll.openmc_mesh_filter_set_translation.restype = c_int
_dll.openmc_mesh_filter_set_translation.errcheck = _error_handler
_dll.openmc_meshborn_filter_get_mesh.argtypes = [c_int32, POINTER(c_int32)]
_dll.openmc_meshborn_filter_get_mesh.restype = c_int
_dll.openmc_meshborn_filter_get_mesh.errcheck = _error_handler
_dll.openmc_meshborn_filter_set_mesh.argtypes = [c_int32, c_int32]
_dll.openmc_meshborn_filter_set_mesh.restype = c_int
_dll.openmc_meshborn_filter_set_mesh.errcheck = _error_handler
_dll.openmc_meshborn_filter_get_translation.argtypes = [c_int32, POINTER(c_double*3)]
_dll.openmc_meshborn_filter_get_translation.restype = c_int
_dll.openmc_meshborn_filter_get_translation.errcheck = _error_handler
_dll.openmc_meshborn_filter_set_translation.argtypes = [c_int32, POINTER(c_double*3)]
_dll.openmc_meshborn_filter_set_translation.restype = c_int
_dll.openmc_meshborn_filter_set_translation.errcheck = _error_handler
_dll.openmc_meshsurface_filter_get_mesh.argtypes = [c_int32, POINTER(c_int32)]
_dll.openmc_meshsurface_filter_get_mesh.restype = c_int
_dll.openmc_meshsurface_filter_get_mesh.errcheck = _error_handler
_dll.openmc_meshsurface_filter_set_mesh.argtypes = [c_int32, c_int32]
_dll.openmc_meshsurface_filter_set_mesh.restype = c_int
_dll.openmc_meshsurface_filter_set_mesh.errcheck = _error_handler
_dll.openmc_meshsurface_filter_get_translation.argtypes = [c_int32, POINTER(c_double*3)]
_dll.openmc_meshsurface_filter_get_translation.restype = c_int
_dll.openmc_meshsurface_filter_get_translation.errcheck = _error_handler
_dll.openmc_meshsurface_filter_set_translation.argtypes = [c_int32, POINTER(c_double*3)]
_dll.openmc_meshsurface_filter_set_translation.restype = c_int
_dll.openmc_meshsurface_filter_set_translation.errcheck = _error_handler
_dll.openmc_new_filter.argtypes = [c_char_p, POINTER(c_int32)]
_dll.openmc_new_filter.restype = c_int
_dll.openmc_new_filter.errcheck = _error_handler
_dll.openmc_spatial_legendre_filter_get_order.argtypes = [c_int32, POINTER(c_int)]
_dll.openmc_spatial_legendre_filter_get_order.restype = c_int
_dll.openmc_spatial_legendre_filter_get_order.errcheck = _error_handler
_dll.openmc_spatial_legendre_filter_set_order.argtypes = [c_int32, c_int]
_dll.openmc_spatial_legendre_filter_set_order.restype = c_int
_dll.openmc_spatial_legendre_filter_set_order.errcheck = _error_handler
_dll.openmc_sphharm_filter_get_order.argtypes = [c_int32, POINTER(c_int)]
_dll.openmc_sphharm_filter_get_order.restype = c_int
_dll.openmc_sphharm_filter_get_order.errcheck = _error_handler
_dll.openmc_sphharm_filter_set_order.argtypes = [c_int32, c_int]
_dll.openmc_sphharm_filter_set_order.restype = c_int
_dll.openmc_sphharm_filter_set_order.errcheck = _error_handler
_dll.openmc_zernike_filter_get_order.argtypes = [c_int32, POINTER(c_int)]
_dll.openmc_zernike_filter_get_order.restype = c_int
_dll.openmc_zernike_filter_get_order.errcheck = _error_handler
_dll.openmc_zernike_filter_set_order.argtypes = [c_int32, c_int]
_dll.openmc_zernike_filter_set_order.restype = c_int
_dll.openmc_zernike_filter_set_order.errcheck = _error_handler
_dll.tally_filters_size.restype = c_size_t

class Filter(_FortranObjectWithID):
    __instances = WeakValueDictionary()

    def __new__(cls, obj=None, uid=None, new=True, index=None):
        mapping = filters
        if index is None:
            if new:
                # Determine ID to assign
                if uid is None:
                    uid = max(mapping, default=0) + 1
                else:
                    if uid in mapping:
                        raise AllocationError('A filter with ID={} has already '
                                              'been allocated.'.format(uid))

                # Set the filter type -- note that the filter_type attribute
                # only exists on subclasses!
                index = c_int32()
                _dll.openmc_new_filter(cls.filter_type.encode(), index)
                index = index.value
            else:
                index = mapping[uid]._index

        if index not in cls.__instances:
            instance = super().__new__(cls)
            instance._index = index
            if uid is not None:
                instance.id = uid
            cls.__instances[index] = instance

        return cls.__instances[index]

    @property
    def id(self):
        filter_id = c_int32()
        _dll.openmc_filter_get_id(self._index, filter_id)
        return filter_id.value

    @id.setter
    def id(self, filter_id):
        _dll.openmc_filter_set_id(self._index, filter_id)

    @property
    def n_bins(self):
        n = c_int()
        _dll.openmc_filter_get_num_bins(self._index, n)
        return n.value


class EnergyFilter(Filter):
    filter_type = 'energy'

    def __init__(self, bins=None, uid=None, new=True, index=None):
        super().__init__(uid, new, index)
        if bins is not None:
            self.bins = bins

    @property
    def bins(self):
        energies = POINTER(c_double)()
        n = c_size_t()
        _dll.openmc_energy_filter_get_bins(self._index, energies, n)
        return as_array(energies, (n.value,))

    @bins.setter
    def bins(self, bins):
        # Get numpy array as a double*
        energies = np.asarray(bins)
        energies_p = energies.ctypes.data_as(POINTER(c_double))

        _dll.openmc_energy_filter_set_bins(
            self._index, len(energies), energies_p)


class CollisionFilter(Filter):
    filter_type = 'collision'


class EnergyoutFilter(EnergyFilter):
    filter_type = 'energyout'


class AzimuthalFilter(Filter):
    filter_type = 'azimuthal'


class CellFilter(Filter):
    filter_type = 'cell'

    @property
    def bins(self):
        cells = POINTER(c_int32)()
        n = c_int32()
        _dll.openmc_cell_filter_get_bins(self._index, cells, n)
        return as_array(cells, (n.value,))


class CellbornFilter(Filter):
    filter_type = 'cellborn'


class CellfromFilter(Filter):
    filter_type = 'cellfrom'


class CellInstanceFilter(Filter):
    filter_type = 'cellinstance'


class DelayedGroupFilter(Filter):
    filter_type = 'delayedgroup'


class DistribcellFilter(Filter):
    filter_type = 'distribcell'


class EnergyFunctionFilter(Filter):
    filter_type = 'energyfunction'

    def __new__(cls, energy=None, y=None, uid=None, new=True, index=None):
        return super().__new__(cls, uid=uid, new=new, index=index)

    def __init__(self, energy=None, y=None, uid=None, new=True, index=None):
        if (energy is None) != (y is None):
            raise AttributeError("Need both energy and y or neither")
        super().__init__(uid, new, index)
        if energy is not None:
            self.set_data(energy, y)

    def set_data(self, energy, y):
        """Set the interpolation information for the filter

        Parameters
        ----------
        energy : numpy.ndarray
            Independent variable for the interpolation
        y : numpy.ndarray
            Dependent variable for the interpolation
        interpolation : {'histogram', 'linear-linear', 'linear-log', 'log-linear', 'log-log', 'quadratic', 'cubic'}
            Interpolation scheme
        """
        energy_array = np.asarray(energy)
        y_array = np.asarray(y)
        energy_p = energy_array.ctypes.data_as(POINTER(c_double))
        y_p = y_array.ctypes.data_as(POINTER(c_double))

        _dll.openmc_energyfunc_filter_set_data(
            self._index, len(energy_array), energy_p, y_p)

    @property
    def energy(self):
        return self._get_attr(_dll.openmc_energyfunc_filter_get_energy)

    @property
    def y(self):
        return self._get_attr(_dll.openmc_energyfunc_filter_get_y)

    @property
    def interpolation(self) -> str:
        interp = c_int()
        _dll.openmc_energyfunc_filter_get_interpolation(self._index, interp)
        return INTERPOLATION_SCHEME[interp.value]

    @interpolation.setter
    def interpolation(self, interp: str):
        interp_ptr = c_char_p(interp.encode())
        _dll.openmc_energyfunc_filter_set_interpolation(self._index, interp_ptr)

    def _get_attr(self, cfunc):
        array_p = POINTER(c_double)()
        n = c_size_t()
        cfunc(self._index, n, array_p)
        return as_array(array_p, (n.value, ))


class LegendreFilter(Filter):
    filter_type = 'legendre'

    def __init__(self, order=None, uid=None, new=True, index=None):
        super().__init__(uid, new, index)
        if order is not None:
            self.order = order

    @property
    def order(self):
        temp_order = c_int()
        _dll.openmc_legendre_filter_get_order(self._index, temp_order)
        return temp_order.value

    @order.setter
    def order(self, order):
        _dll.openmc_legendre_filter_set_order(self._index, order)


class MaterialFilter(Filter):
    filter_type = 'material'

    def __init__(self, bins=None, uid=None, new=True, index=None):
        super().__init__(uid, new, index)
        if bins is not None:
            self.bins = bins

    @property
    def bins(self):
        materials = POINTER(c_int32)()
        n = c_size_t()
        _dll.openmc_material_filter_get_bins(self._index, materials, n)
        return [Material(index=materials[i]) for i in range(n.value)]

    @bins.setter
    def bins(self, materials):
        # Get material indices as int32_t[]
        n = len(materials)
        bins = (c_int32*n)(*(m._index for m in materials))
        _dll.openmc_material_filter_set_bins(self._index, n, bins)


class MaterialFromFilter(Filter):
    filter_type = 'materialfrom'


class MeshFilter(Filter):
    """Mesh filter stored internally.

    This class exposes a Mesh filter that is stored internally in the OpenMC
    library. To obtain a view of a Mesh filter with a given ID, use the
    :data:`openmc.lib.filters` mapping.

    Parameters
    ----------
    mesh : openmc.lib.Mesh
        Mesh to use for the filter
    uid : int or None
        Unique ID of the Mesh filter
    new : bool
        When `index` is None, this argument controls whether a new object is
        created or a view of an existing object is returned.
    index : int
        Index in the `filters` array.

    Attributes
    ----------
    filter_type : str
        Type of filter
    mesh : openmc.lib.Mesh
        Mesh used for the filter
    translation : Iterable of float
        3-D coordinates of the translation vector

    """
    filter_type = 'mesh'

    def __init__(self, mesh=None, uid=None, new=True, index=None):
        super().__init__(uid, new, index)
        if mesh is not None:
            self.mesh = mesh

    @property
    def mesh(self):
        index_mesh = c_int32()
        _dll.openmc_mesh_filter_get_mesh(self._index, index_mesh)
        return _get_mesh(index_mesh.value)

    @mesh.setter
    def mesh(self, mesh):
        _dll.openmc_mesh_filter_set_mesh(self._index, mesh._index)

    @property
    def translation(self):
        translation = (c_double*3)()
        _dll.openmc_mesh_filter_get_translation(self._index, translation)
        return tuple(translation)

    @translation.setter
    def translation(self, translation):
        _dll.openmc_mesh_filter_set_translation(self._index, (c_double*3)(*translation))


class MeshBornFilter(Filter):
    """MeshBorn filter stored internally.

    This class exposes a MeshBorn filter that is stored internally in the OpenMC
    library. To obtain a view of a MeshBorn filter with a given ID, use the
    :data:`openmc.lib.filters` mapping.

    Parameters
    ----------
    mesh : openmc.lib.Mesh
        Mesh to use for the filter
    uid : int or None
        Unique ID of the MeshBorn filter
    new : bool
        When `index` is None, this argument controls whether a new object is
        created or a view of an existing object is returned.
    index : int
        Index in the `filters` array.

    Attributes
    ----------
    filter_type : str
        Type of filter
    mesh : openmc.lib.Mesh
        Mesh used for the filter
    translation : Iterable of float
        3-D coordinates of the translation vector

    """
    filter_type = 'meshborn'

    def __init__(self, mesh=None, uid=None, new=True, index=None):
        super().__init__(uid, new, index)
        if mesh is not None:
            self.mesh = mesh

    @property
    def mesh(self):
        index_mesh = c_int32()
        _dll.openmc_meshborn_filter_get_mesh(self._index, index_mesh)
        return _get_mesh(index_mesh.value)

    @mesh.setter
    def mesh(self, mesh):
        _dll.openmc_meshborn_filter_set_mesh(self._index, mesh._index)

    @property
    def translation(self):
        translation = (c_double*3)()
        _dll.openmc_meshborn_filter_get_translation(self._index, translation)
        return tuple(translation)

    @translation.setter
    def translation(self, translation):
        _dll.openmc_meshborn_filter_set_translation(self._index, (c_double*3)(*translation))


class MeshSurfaceFilter(Filter):
    """MeshSurface filter stored internally.

    This class exposes a MeshSurface filter that is stored internally in the
    OpenMC library. To obtain a view of a MeshSurface filter with a given ID,
    use the :data:`openmc.lib.filters` mapping.

    Parameters
    ----------
    mesh : openmc.lib.Mesh
        Mesh to use for the filter
    uid : int or None
        Unique ID of the MeshSurface filter
    new : bool
        When `index` is None, this argument controls whether a new object is
        created or a view of an existing object is returned.
    index : int
        Index in the `filters` array.

    Attributes
    ----------
    filter_type : str
        Type of filter
    mesh : openmc.lib.Mesh
        Mesh used for the filter
    translation : Iterable of float
        3-D coordinates of the translation vector

    """
    filter_type = 'meshsurface'

    def __init__(self, mesh=None, uid=None, new=True, index=None):
        super().__init__(uid, new, index)
        if mesh is not None:
            self.mesh = mesh

    @property
    def mesh(self):
        index_mesh = c_int32()
        _dll.openmc_meshsurface_filter_get_mesh(self._index, index_mesh)
        return _get_mesh(index_mesh.value)

    @mesh.setter
    def mesh(self, mesh):
        _dll.openmc_meshsurface_filter_set_mesh(self._index, mesh._index)

    @property
    def translation(self):
        translation = (c_double*3)()
        _dll.openmc_meshsurface_filter_get_translation(self._index, translation)
        return tuple(translation)

    @translation.setter
    def translation(self, translation):
        _dll.openmc_meshsurface_filter_set_translation(self._index, (c_double*3)(*translation))


class MuFilter(Filter):
    filter_type = 'mu'


class ParticleFilter(Filter):
    filter_type = 'particle'

    @property
    def bins(self):
        particle_i = np.zeros((self.n_bins,), dtype=c_int)
        _dll.openmc_particle_filter_get_bins(
            self._index, particle_i.ctypes.data_as(POINTER(c_int)))
        return [ParticleType(i) for i in particle_i]


class PolarFilter(Filter):
    filter_type = 'polar'


class SphericalHarmonicsFilter(Filter):
    filter_type = 'sphericalharmonics'

    def __init__(self, order=None, uid=None, new=True, index=None):
        super().__init__(uid, new, index)
        if order is not None:
            self.order = order

    @property
    def order(self):
        temp_order = c_int()
        _dll.openmc_sphharm_filter_get_order(self._index, temp_order)
        return temp_order.value

    @order.setter
    def order(self, order):
        _dll.openmc_sphharm_filter_set_order(self._index, order)


class SpatialLegendreFilter(Filter):
    filter_type = 'spatiallegendre'

    def __init__(self, order=None, uid=None, new=True, index=None):
        super().__init__(uid, new, index)
        if order is not None:
            self.order = order

    @property
    def order(self):
        temp_order = c_int()
        _dll.openmc_spatial_legendre_filter_get_order(self._index, temp_order)
        return temp_order.value

    @order.setter
    def order(self, order):
        _dll.openmc_spatial_legendre_filter_set_order(self._index, order)


class SurfaceFilter(Filter):
    filter_type = 'surface'


class UniverseFilter(Filter):
    filter_type = 'universe'


class ZernikeFilter(Filter):
    filter_type = 'zernike'

    def __init__(self, order=None, uid=None, new=True, index=None):
        super().__init__(uid, new, index)
        if order is not None:
            self.order = order

    @property
    def order(self):
        temp_order = c_int()
        _dll.openmc_zernike_filter_get_order(self._index, temp_order)
        return temp_order.value

    @order.setter
    def order(self, order):
        _dll.openmc_zernike_filter_set_order(self._index, order)


class ZernikeRadialFilter(ZernikeFilter):
    filter_type = 'zernikeradial'


_FILTER_TYPE_MAP = {
    'azimuthal': AzimuthalFilter,
    'cell': CellFilter,
    'cellborn': CellbornFilter,
    'cellfrom': CellfromFilter,
    'cellinstance': CellInstanceFilter,
    'delayedgroup': DelayedGroupFilter,
    'distribcell': DistribcellFilter,
    'energy': EnergyFilter,
    'energyout': EnergyoutFilter,
    'energyfunction': EnergyFunctionFilter,
    'legendre': LegendreFilter,
    'material': MaterialFilter,
    'materialfrom': MaterialFromFilter,
    'mesh': MeshFilter,
    'meshborn': MeshBornFilter,
    'meshsurface': MeshSurfaceFilter,
    'mu': MuFilter,
    'particle': ParticleFilter,
    'polar': PolarFilter,
    'sphericalharmonics': SphericalHarmonicsFilter,
    'spatiallegendre': SpatialLegendreFilter,
    'surface': SurfaceFilter,
    'time': TimeFilter,
    'timedmesh': TimedMeshFilter,
    'universe': UniverseFilter,
    'zernike': ZernikeFilter,
    'zernikeradial': ZernikeRadialFilter
}


def _get_filter(index):
    filter_type = create_string_buffer(20)
    _dll.openmc_filter_get_type(index, filter_type)
    filter_type = filter_type.value.decode()
    return _FILTER_TYPE_MAP[filter_type](index=index)


class _FilterMapping(Mapping):
    def __getitem__(self, key):
        index = c_int32()
        try:
            _dll.openmc_get_filter_index(key, index)
        except (AllocationError, InvalidIDError) as e:
            # __contains__ expects a KeyError to work correctly
            raise KeyError(str(e))
        return _get_filter(index.value)

    def __iter__(self):
        for i in range(len(self)):
            yield _get_filter(i).id

    def __len__(self):
        return _dll.tally_filters_size()

    def __repr__(self):
        return repr(dict(self))

filters = _FilterMapping()
