from collections import Mapping
from ctypes import c_int, c_int32, c_double, c_char_p, POINTER, \
    create_string_buffer
from weakref import WeakValueDictionary

import numpy as np
from numpy.ctypeslib import as_array

from . import _dll
from .core import _FortranObjectWithID
from .error import _error_handler, AllocationError, InvalidIDError
from .material import Material


__all__ = ['Filter', 'AzimuthalFilter', 'CellFilter',
           'CellbornFilter', 'CellfromFilter', 'DistribcellFilter',
           'DelayedGroupFilter', 'EnergyFilter', 'EnergyoutFilter',
           'EnergyFunctionFilter', 'MaterialFilter', 'MeshFilter',
           'MuFilter', 'PolarFilter', 'SurfaceFilter',
           'UniverseFilter', 'filters']

# Tally functions
_dll.openmc_energy_filter_get_bins.argtypes = [
    c_int32, POINTER(POINTER(c_double)), POINTER(c_int32)]
_dll.openmc_energy_filter_get_bins.restype = c_int
_dll.openmc_energy_filter_get_bins.errcheck = _error_handler
_dll.openmc_energy_filter_set_bins.argtypes = [c_int32, c_int32, POINTER(c_double)]
_dll.openmc_energy_filter_set_bins.restype = c_int
_dll.openmc_energy_filter_set_bins.errcheck = _error_handler
_dll.openmc_extend_filters.argtypes = [c_int32, POINTER(c_int32), POINTER(c_int32)]
_dll.openmc_extend_filters.restype = c_int
_dll.openmc_extend_filters.errcheck = _error_handler
_dll.openmc_filter_get_id.argtypes = [c_int32, POINTER(c_int32)]
_dll.openmc_filter_get_id.restype = c_int
_dll.openmc_filter_get_id.errcheck = _error_handler
_dll.openmc_filter_get_type.argtypes = [c_int32, c_char_p]
_dll.openmc_filter_get_type.restype = c_int
_dll.openmc_filter_get_type.errcheck = _error_handler
_dll.openmc_filter_set_id.argtypes = [c_int32, c_int32]
_dll.openmc_filter_set_id.restype = c_int
_dll.openmc_filter_set_id.errcheck = _error_handler
_dll.openmc_filter_set_type.argtypes = [c_int32, c_char_p]
_dll.openmc_filter_set_type.restype = c_int
_dll.openmc_filter_set_type.errcheck = _error_handler
_dll.openmc_get_filter_index.argtypes = [c_int32, POINTER(c_int32)]
_dll.openmc_get_filter_index.restype = c_int
_dll.openmc_get_filter_index.errcheck = _error_handler
_dll.openmc_material_filter_get_bins.argtypes = [
    c_int32, POINTER(POINTER(c_int32)), POINTER(c_int32)]
_dll.openmc_material_filter_get_bins.restype = c_int
_dll.openmc_material_filter_get_bins.errcheck = _error_handler
_dll.openmc_material_filter_set_bins.argtypes = [c_int32, c_int32, POINTER(c_int32)]
_dll.openmc_material_filter_set_bins.restype = c_int
_dll.openmc_material_filter_set_bins.errcheck = _error_handler
_dll.openmc_mesh_filter_set_mesh.argtypes = [c_int32, c_int32]
_dll.openmc_mesh_filter_set_mesh.restype = c_int
_dll.openmc_mesh_filter_set_mesh.errcheck = _error_handler


class Filter(_FortranObjectWithID):
    __instances = WeakValueDictionary()

    def __new__(cls, obj=None, uid=None, new=True, index=None):
        mapping = filters
        if index is None:
            if new:
                # Determine ID to assign
                if uid is None:
                    try:
                        uid = max(mapping) + 1
                    except ValueError:
                        uid = 1
                else:
                    if uid in mapping:
                        raise AllocationError('A filter with ID={} has already '
                                              'been allocated.'.format(uid))

                # Resize internal array
                index = c_int32()
                _dll.openmc_extend_filters(1, index, None)

                # Set the filter type -- note that the filter_type attribute
                # only exists on subclasses!
                _dll.openmc_filter_set_type(index, cls.filter_type.encode())
                index = index.value
            else:
                index = mapping[uid]._index

        if index not in cls.__instances:
            instance = super(Filter, cls).__new__(cls)
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


class EnergyFilter(Filter):
    filter_type = 'energy'

    def __init__(self, bins=None, uid=None, new=True, index=None):
        super(EnergyFilter, self).__init__(uid, new, index)
        if bins is not None:
            self.bins = bins

    @property
    def bins(self):
        energies = POINTER(c_double)()
        n = c_int32()
        _dll.openmc_energy_filter_get_bins(self._index, energies, n)
        return as_array(energies, (n.value,))

    @bins.setter
    def bins(self, bins):
        # Get numpy array as a double*
        energies = np.asarray(bins)
        energies_p = energies.ctypes.data_as(POINTER(c_double))

        _dll.openmc_energy_filter_set_bins(
            self._index, len(energies), energies_p)


class EnergyoutFilter(Filter):
    filter_type = 'energyout'


class AzimuthalFilter(Filter):
    filter_type = 'azimuthal'


class CellFilter(Filter):
    filter_type = 'cell'


class CellbornFilter(Filter):
    filter_type = 'cellborn'


class CellfromFilter(Filter):
    filter_type = 'cellfrom'


class DelayedGroupFilter(Filter):
    filter_type = 'delayedgroup'


class DistribcellFilter(Filter):
    filter_type = 'distribcell'


class EnergyFunctionFilter(Filter):
    filter_type = 'energyfunction'


class MaterialFilter(Filter):
    filter_type = 'material'

    def __init__(self, bins=None, uid=None, new=True, index=None):
        super(MaterialFilter, self).__init__(uid, new, index)
        if bins is not None:
            self.bins = bins

    @property
    def bins(self):
        materials = POINTER(c_int32)()
        n = c_int32()
        _dll.openmc_material_filter_get_bins(self._index, materials, n)
        return [Material(index=materials[i]) for i in range(n.value)]

    @bins.setter
    def bins(self, materials):
        # Get material indices as int32_t[]
        n = len(materials)
        bins = (c_int32*n)(*(m._index for m in materials))

        _dll.openmc_material_filter_set_bins(self._index, n, bins)


class MeshFilter(Filter):
    filter_type = 'mesh'


class MuFilter(Filter):
    filter_type = 'mu'


class PolarFilter(Filter):
    filter_type = 'polar'


class SurfaceFilter(Filter):
    filter_type = 'surface'


class UniverseFilter(Filter):
    filter_type = 'universe'


_FILTER_TYPE_MAP = {
    'azimuthal': AzimuthalFilter,
    'cell': CellFilter,
    'cellborn': CellbornFilter,
    'cellfrom': CellfromFilter,
    'delayedgroup': DelayedGroupFilter,
    'distribcell': DistribcellFilter,
    'energy': EnergyFilter,
    'energyout': EnergyoutFilter,
    'energyfunction': EnergyFunctionFilter,
    'material': MaterialFilter,
    'mesh': MeshFilter,
    'mu': MuFilter,
    'polar': PolarFilter,
    'surface': SurfaceFilter,
    'universe': UniverseFilter,
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
            yield _get_filter(i + 1).id

    def __len__(self):
        return c_int32.in_dll(_dll, 'n_filters').value

    def __repr__(self):
        return repr(dict(self))

filters = _FilterMapping()
