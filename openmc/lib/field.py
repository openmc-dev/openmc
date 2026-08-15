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
_dll.openmc_temperature_field_size.restype = c_size_t


class TemperatureField:
    """Temperature field stored internally in the OpenMC library."""

    def __len__(self):
        return self.size

    @property
    def size(self):
        """Number of cells in the temperature field."""
        return _dll.openmc_temperature_field_size()

    def __getitem__(self, index):
        temperature = c_double()
        _dll.openmc_temperature_field_get_value(index, temperature)
        return temperature.value

    def __setitem__(self, index, temperature):
        _dll.openmc_temperature_field_set_temperature(index, temperature)


temperature_field = TemperatureField()
