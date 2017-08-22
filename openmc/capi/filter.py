from collections import Mapping
from ctypes import c_int, c_int32, c_double, c_char_p, POINTER, \
    create_string_buffer
from weakref import WeakValueDictionary

import numpy as np
from numpy.ctypeslib import as_array

from . import _dll
from .error import _error_handler
from .material import MaterialView


__all__ = ['FilterView', 'AzimuthalFilterView', 'CellFilterView',
           'CellbornFilterView', 'CellfromFilterView', 'DistribcellFilterView',
           'DelayedGroupFilterView', 'EnergyFilterView', 'EnergyoutFilterView',
           'EnergyFunctionFilterView', 'MaterialFilterView', 'MeshFilterView',
           'MuFilterView', 'PolarFilterView', 'SurfaceFilterView',
           'UniverseFilterView', 'filters']

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
_dll.openmc_get_filter.argtypes = [c_int32, POINTER(c_int32)]
_dll.openmc_get_filter.restype = c_int
_dll.openmc_get_filter.errcheck = _error_handler
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


class FilterView(object):
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
        filter_id = c_int32()
        _dll.openmc_filter_get_id(self._index, filter_id)
        return filter_id.value

    @id.setter
    def id(self, filter_id):
        _dll.openmc_filter_set_id(self._index, filter_id)


class EnergyFilterView(FilterView):
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


class EnergyoutFilterView(FilterView):
    pass


class AzimuthalFilterView(FilterView):
    pass


class CellFilterView(FilterView):
    pass


class CellbornFilterView(FilterView):
    pass


class CellfromFilterView(FilterView):
    pass


class DelayedGroupFilterView(FilterView):
    pass


class DistribcellFilterView(FilterView):
    pass


class EnergyFunctionFilterView(FilterView):
    pass


class MaterialFilterView(FilterView):
    @property
    def bins(self):
        materials = POINTER(c_int32)()
        n = c_int32()
        _dll.openmc_material_filter_get_bins(self._index, materials, n)
        return [MaterialView(materials[i]) for i in range(n.value)]

    @bins.setter
    def bins(self, materials):
        # Get material indices as int32_t[]
        n = len(materials)
        bins = (c_int32*n)(*(m._index for m in materials))

        _dll.openmc_material_filter_set_bins(self._index, n, bins)

    @classmethod
    def new(cls, bins=None):
        index = c_int32()
        _dll.openmc_extend_filters(1, index, None)
        _dll.openmc_filter_set_type(index, b'material')
        f = cls(index.value)
        if bins is not None:
            f.bins = bins
        return f


class MeshFilterView(FilterView):
    pass


class MuFilterView(FilterView):
    pass


class PolarFilterView(FilterView):
    pass


class SurfaceFilterView(FilterView):
    pass


class UniverseFilterView(FilterView):
    pass


_filter_type_map = {
    'azimuthal': AzimuthalFilterView,
    'cell': CellFilterView,
    'cellborn': CellbornFilterView,
    'cellfrom': CellfromFilterView,
    'delayedgroup': DelayedGroupFilterView,
    'distribcell': DistribcellFilterView,
    'energy': EnergyFilterView,
    'energyout': EnergyoutFilterView,
    'energyfunction': EnergyFunctionFilterView,
    'material': MaterialFilterView,
    'mesh': MeshFilterView,
    'mu': MuFilterView,
    'polar': PolarFilterView,
    'surface': SurfaceFilterView,
    'universe': UniverseFilterView,
}


def _get_filter(index):
    filter_type = create_string_buffer(20)
    _dll.openmc_filter_get_type(index, filter_type)
    filter_type = filter_type.value.decode()
    return _filter_type_map[filter_type](index)


class _FilterMapping(Mapping):
    def __getitem__(self, key):
        index = c_int32()
        _dll.openmc_get_filter(key, index)
        return _get_filter(index.value)

    def __iter__(self):
        for i in range(len(self)):
            yield TallyView(i + 1).id

    def __len__(self):
        return c_int32.in_dll(_dll, 'n_filters').value

filters = _FilterMapping()
