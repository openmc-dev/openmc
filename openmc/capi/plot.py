from ctypes import c_int, c_int32, c_double, c_ushort, POINTER, Structure

from . import _dll
from .core import _DLLGlobal

import numpy as np

class _Position(Structure):
    pass

_Position._fields_ = [('x', c_double),
                      ('y', c_double),
                      ('z', c_double)]

class _RGBColor(Structure):
    pass

_RGBColor._fields_ = [('red', c_ushort),
                      ('green', c_ushort),
                      ('blue', c_ushort)]

class _Plot(Structure):
    pass

_Plot._fields_ = [('origin_', _Position),
                  ('width_', _Position),
                  ('basis_', c_int),
                  ('pixels_', c_int*3),
                  ('level_', c_int)]


_dll.openmc_id_map.argtypes= [POINTER(_Plot),]
_dll.openmc_id_map.restype = c_int

def image_data_for_plot(plot):
    img_data = np.zeros((plot.pixels_[0], plot.pixels_[1], 2), dtype=np.dtype('int32'))
    out = _dll.openmc_id_map(POINTER(_Plot)(plot), img_data.ctypes.data_as(POINTER(c_int32)))
    return img_data
