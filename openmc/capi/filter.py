from ctypes import c_int, c_int32, c_double, c_char_p, POINTER

from . import _dll
from .error import _error_handler


__all__ = []

# Tally functions
_dll.openmc_extend_filters.argtypes = [c_int32, POINTER(c_int32), POINTER(c_int32)]
_dll.openmc_extend_filters.restype = c_int
_dll.openmc_extend_filters.errcheck = _error_handler
_dll.openmc_filter_set_type.argtypes = [c_int32, c_char_p]
_dll.openmc_filter_set_type.restype = c_int
_dll.openmc_filter_set_type.errcheck = _error_handler
