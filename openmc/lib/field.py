from collections.abc import Sequence
from ctypes import c_double, c_int32, POINTER
import operator

from . import _dll
from .error import _error_handler
from .mesh import _get_mesh

__all__ = ['TemperatureField', 'temperature_field']


_dll.openmc_temperature_field_get_temperature.argtypes = [
    c_int32, POINTER(c_double)]
_dll.openmc_temperature_field_get_temperature.restype = c_int32
_dll.openmc_temperature_field_get_temperature.errcheck = _error_handler
_dll.openmc_temperature_field_set_temperature.argtypes = [c_int32, c_double]
_dll.openmc_temperature_field_set_temperature.restype = c_int32
_dll.openmc_temperature_field_set_temperature.errcheck = _error_handler
_dll.openmc_temperature_field_get_mesh.argtypes = [POINTER(c_int32)]
_dll.openmc_temperature_field_get_mesh.restype = c_int32
_dll.openmc_temperature_field_get_mesh.errcheck = _error_handler


class TemperatureField(Sequence):
    """Mutable sequence of temperatures stored in the OpenMC library.

    Element *i* holds the temperature in [K] of element *i* of the mesh that
    the temperature field is defined on. Temperatures may be changed between
    batches; the new value must lie within the temperature range covered by
    the loaded cross sections.

    """

    def __len__(self):
        mesh = self.mesh
        return 0 if mesh is None else mesh.n_elements

    @property
    def mesh(self):
        """Mesh associated with the temperature field, or None."""
        index = c_int32()
        _dll.openmc_temperature_field_get_mesh(index)
        return None if index.value < 0 else _get_mesh(index.value)

    def _index(self, index):
        """Normalize an index, converting Python conventions to C ones."""
        index = operator.index(index)
        size = len(self)
        if index < 0:
            index += size
        if index < 0 or index >= size:
            raise IndexError('Index in temperature field is out of bounds.')
        return index

    def __getitem__(self, index):
        if isinstance(index, slice):
            return [self[i] for i in range(*index.indices(len(self)))]

        temperature = c_double()
        _dll.openmc_temperature_field_get_temperature(
            self._index(index), temperature)
        return temperature.value

    def __setitem__(self, index, temperature):
        _dll.openmc_temperature_field_set_temperature(
            self._index(index), temperature)


temperature_field = TemperatureField()
