from collections.abc import Sequence
from ctypes import c_double, c_int32, c_size_t, POINTER

from . import _dll
from .error import _error_handler

__all__ = ['TemperatureField', 'temperature_field']


_dll.openmc_temperature_field_get_value.argtypes = [
    c_int32, POINTER(c_double)]
_dll.openmc_temperature_field_get_value.restype = c_int32
_dll.openmc_temperature_field_get_value.errcheck = _error_handler
_dll.openmc_temperature_field_set_temperature.argtypes = [c_int32, c_double]
_dll.openmc_temperature_field_set_temperature.restype = c_int32
_dll.openmc_temperature_field_set_temperature.errcheck = _error_handler
_dll.openmc_temperature_field_size.argtypes = []
_dll.openmc_temperature_field_size.restype = c_size_t


class TemperatureField(Sequence):
    """Mutable sequence of temperatures stored in the OpenMC library.

    Element *i* holds the temperature in [K] of element *i* of the mesh that
    the temperature field is defined on. Temperatures may be changed between
    batches; the new value must lie within the temperature range covered by
    the loaded cross sections.

    """

    def __len__(self):
        return self.size

    @property
    def size(self):
        """Number of elements in the temperature field."""
        return _dll.openmc_temperature_field_size()

    def _index(self, index):
        """Normalize an index, converting Python conventions to C ones."""
        index = index.__index__()
        if index < 0:
            index += self.size
        if index < 0 or index >= self.size:
            raise IndexError('Index in temperature field is out of bounds.')
        return index

    def __getitem__(self, index):
        temperature = c_double()
        _dll.openmc_temperature_field_get_value(
            self._index(index), temperature)
        return temperature.value

    def __setitem__(self, index, temperature):
        _dll.openmc_temperature_field_set_temperature(
            self._index(index), temperature)


temperature_field = TemperatureField()
