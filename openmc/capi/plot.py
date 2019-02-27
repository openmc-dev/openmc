from ctypes import c_int, c_int32, c_double, c_ushort, POINTER, Structure

from . import _dll
from .core import _DLLGlobal
from .error import _error_handler

import numpy as np

class _IntThree():
    type_ = c_int*3
    val_ = None
    def __init__(self, vals=None):
        if vals and iterable(vals) and all(vals == int):
            self.val_ = type_(vals[0], vals[1], vals[2])
        else:
            raise ValueError("{} is not a valid object to construct IntThree.".format(vals))

    def __str__(self):
        return "({}, {}, {})".format(self.val_[0], self.val_[1], self.val_[2])

class _Position(Structure):
    _fields_ = [('x', c_double),
                ('y', c_double),
                ('z', c_double)]

    def __init__(self, vals=None):
        if vals:
            x = vals[0]
            y = vals[1]
            z = vals[2]
        else:
            x = 0.0
            y = 0.0
            z = 0.0

    @property
    def x(self):
        return self.x

    @property
    def y(self):
        return self.y

    @property
    def z(self):
        return self.z

    @x.setter
    def x(self, x_val):
        assert(isinstance(x_val, float))
        self.x = x_val

    @y.setter
    def y(self, y_val):
        assert(isinstance(y_val, float))
        self.y = y_val

    @z.setter
    def z(self, z_val):
        assert(isinstance(z_val, float))
        self.z = z_val

    def __str__(self):
        return "({}, {}, {})".format(self.x, self.y, self.z)

class _Plot(Structure):
    _fields_ = [('origin_', _Position),
                ('width_', _Position),
                ('basis_', c_int),
                ('pixels_', c_int*3),
                ('level_', c_int)]

    def __init__(self):
        self.level_ = -1

    @property
    def origin(self):
        out = [self.origin_.x, self.origin_.y, self.origin_.z]
        return [float(val) for val in out]

    @property
    def width(self):
        return float(self.width_.x)

    @property
    def height(self):
        return float(self.width_.y)

    @property
    def basis(self):
        if self.basis_ == 1:
            return 'xy'
        elif self.basis_ == 2:
            return 'xz'
        elif self.basis_ == 3:
            return 'yz'

        raise ValueError("Plot basis {} is invalid".format(basis_))

    @property
    def hRes(self):
        return self.pixels_[0]

    @property
    def vRes(self):
        return self.pixels_[1]

    @property
    def level(self):
        return int(self.level_)

    @origin.setter
    def origin(self, origin):
        self.origin_.x = origin[0]
        self.origin_.y = origin[1]
        self.origin_.z = origin[2]

    @width.setter
    def width(self, width):
        self.width_.x = width

    @height.setter
    def height(self, height):
        self.width_.y = height

    @basis.setter
    def basis(self, basis):
        if isinstance(basis, str):
            valid_bases = ('xy', 'xz', 'yz')
            if basis not in valid_bases:
                raise ValueError("{} is not a valid plot basis.".format(basis))

            if basis == 'xy':
                self.basis_ = 1
            elif basis == 'xz':
                self.basis_ = 2
            elif basis == 'yz':
                self.basis_ = 3
            return

        if isinstance(basis, int):
            valid_bases = (1, 2, 3)
            if basis not in valid_bases:
                raise ValueError("{} is not a valid plot basis.".format(basis))
            self.basis_ = basis
            return

        raise ValueError("{} of type {} is an invalid plot basis".format(basis, type(basis)))

    @hRes.setter
    def hRes(self, hRes):
        self.pixels_[0] = hRes

    @vRes.setter
    def vRes(self, vRes):
        self.pixels_[1] = vRes

    @level.setter
    def level(self, level):
        self.level_ = level

    def __str__(self):
        out_str = "-------"
        out_str += "Plot:\n"
        out_str += "-------"
        out_str += "Origin: " + str(self.origin) + "\n"
        out_str += "Width: " + str(self.width) + "\n"
        out_str += "Height: " + str(self.height) + "\n"
        out_str += "Basis: " + str(self.basis) + "\n"
        out_str += "HRes: " + str(self.hRes) + "\n"
        out_str += "VRes: " + str(self.vRes) + "\n"
        out_str += "Level: " + str(self.level) + "\n"
        return out_str


_dll.openmc_id_map.argtypes= [POINTER(_Plot),]
_dll.openmc_id_map.restype = c_int
_dll.openmc_id_map.errcheck = _error_handler

def image_data_for_plot(plot):
    img_data = np.zeros((plot.pixels_[0], plot.pixels_[1], 2), dtype=np.dtype('int32'))
    out = _dll.openmc_id_map(POINTER(_Plot)(plot), img_data.ctypes.data_as(POINTER(c_int32)))
    return img_data
