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

    def __new__(cls, filter_type, uid=None, new=True, index=None):
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

                index = c_int32()
                _dll.openmc_extend_filters(1, index, None)
                _dll.openmc_filter_set_type(index, filter_type)
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
    def __new__(cls, bins=None, uid=None, new=True, index=None):
        return super(EnergyFilter, cls).__new__(cls, b'energy', uid, new, index)

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
    def __new__(cls, bins=None, uid=None, new=True, index=None):
        return super(Energyoutfilter, cls).__new__(cls, b'energyout',
                                                   uid, new, index)


class AzimuthalFilter(Filter):
    def __new__(cls, bins=None, uid=None, new=True, index=None):
        return super(AzimuthalFilter, cls).__new__(cls, b'azimuthal',
                                                   uid, new, index)


class CellFilter(Filter):
    def __new__(cls, bins=None, uid=None, new=True, index=None):
        return super(CellFilter, cls).__new__(cls, b'cell', uid, new, index)


class CellbornFilter(Filter):
    def __new__(cls, bins=None, uid=None, new=True, index=None):
        return super(CellbornFilter, cls).__new__(cls, b'cellborn', uid,
                                                  new, index)


class CellfromFilter(Filter):
    def __new__(cls, bins=None, uid=None, new=True, index=None):
        return super(CellfromFilter, cls).__new__(cls, b'cellfrom', uid,
                                                  new, index)


class DelayedGroupFilter(Filter):
    def __new__(cls, bins=None, uid=None, new=True, index=None):
        return super(DelayedGroupFilter, cls).__new__(cls, b'delayedgroup',
                                                      uid, new, index)


class DistribcellFilter(Filter):
    def __new__(cls, bins=None, uid=None, new=True, index=None):
        return super(DistribcellFilter, cls).__new__(cls, b'distribcell',
                                                     uid, new, index)


class EnergyFunctionFilter(Filter):
    def __new__(cls, bins=None, uid=None, new=True, index=None):
        return super(EnergyFunctionFilter, cls).__new__(
            cls, b'energyfunction', uid, new, index)


class MaterialFilter(Filter):
    def __new__(cls, bins=None, uid=None, new=True, index=None):
        return super(MaterialFilter, cls).__new__(cls, b'material',
                                                   uid, new, index)

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
    def __new__(cls, bins=None, uid=None, new=True, index=None):
        return super(MeshFilter, cls).__new__(cls, b'mesh', uid, new, index)


class MuFilter(Filter):
    def __new__(cls, bins=None, uid=None, new=True, index=None):
        return super(MuFilter, cls).__new__(cls, b'mu', uid, new, index)


class PolarFilter(Filter):
    def __new__(cls, bins=None, uid=None, new=True, index=None):
        return super(PolarFilter, cls).__new__(cls, b'polar', uid, new, index)


class SurfaceFilter(Filter):
    def __new__(cls, bins=None, uid=None, new=True, index=None):
        return super(SurfaceFilter, cls).__new__(cls, b'surface', uid,
                                                 new, index)


class UniverseFilter(Filter):
    def __new__(cls, bins=None, uid=None, new=True, index=None):
        return super(UniverseFilter, cls).__new__(cls, b'universe', uid,
                                                  new, index)


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
