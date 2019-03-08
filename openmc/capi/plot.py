from ctypes import c_int, c_int32, c_double, Structure, POINTER

from . import _dll
from .core import _DLLGlobal
from .error import _error_handler

import numpy as np


class _Position(Structure):
    """Definition of an xyz location in space with underlying c-types

    C-type Attributes
    -----------------
    x : c_double
        Position's x value (default: 0.0)
    y : c_double
        Position's y value (default: 0.0)
    z : c_double
        Position's z value (default: 0.0)
    """
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
        return "Position: ({}, {}, {})".format(self.x, self.y, self.z)


class _PlotBase(Structure):
    """A structure defining a 2-D geometry slice with underlying c-types

    C-Type Attributes
    -----------------
    origin : openmc.capi.plot._Position
        A position defining the origin of the plot.
    width_ : openmc.capi.plot._Position
        The width of the plot along the x, y, and z axes, respectively
    basis_ : c_int
        The axes basis of the plot view.
    pixels_ : c_int[3]
        The resolution of the plot in the horizontal and vertical dimensions
    level_ : c_int
        The universe level for the plot view

    Attributes
    ----------
    origin : tuple or list of ndarray
        Origin (center) of the plot
    width : float
        The horizontal dimension of the plot in geometry units (cm)
    height : float
        The vertical dimension of the plot in geometry units (cm)
    basis : string
        One of {'xy', 'xz', 'yz'} indicating the horizontal and vertical
        axes of the plot.
    hRes : float
        The horizontal resolution of the plot in pixels
    vRes : float
        The vertical resolution of the plot in pixels
    level : int
        The universe level for the plot (default: -1 -> all universes shown)
    """
    _fields_ = [('origin_', _Position),
                ('width_', _Position),
                ('basis_', c_int),
                ('pixels_', 3*c_int),
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
            basis = basis.lower()
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

        raise ValueError("{} of type {} is an"
                         " invalid plot basis".format(basis, type(basis)))

    @hRes.setter
    def hRes(self, hRes):
        self.pixels_[0] = hRes

    @vRes.setter
    def vRes(self, vRes):
        self.pixels_[1] = vRes

    @level.setter
    def level(self, level):
        self.level_ = level

    def __repr__(self):
        out_str = "-----\n"
        out_str += "Plot:\n"
        out_str += "-----\n"
        out_str += "Origin: {}\n".format(self.origin)
        out_str += "Width: {}\n".format(self.width)
        out_str += "Height: {}\n".format(self.height)
        out_str += "Basis: {}\n".format(self.basis)
        out_str += "HRes: {}\n".format(self.hRes)
        out_str += "VRes: {}\n".format(self.vRes)
        out_str += "Level: {}\n".format(self.level)
        return out_str

    def __str__(self):
        return self.__repr__()


_dll.openmc_id_map.argtypes = [POINTER(_PlotBase), POINTER(c_int32)]
_dll.openmc_id_map.restype = c_int
_dll.openmc_id_map.errcheck = _error_handler


def id_map(plot):
    """
    Generate a 2-D map of (cell_id, material_id). Used for in-memory image
         generation.

    Parameters
    ----------
    plot : An openmc.capi.plot._PlotBase object describing the slice of the
           model to be generated

    Returns
    -------
    id_map : a NumPy array with shape (vertical pixels, horizontal pixels, 2)
             of OpenMC property ids with dtype int32

    """
    img_data = np.zeros((plot.vRes, plot.hRes, 2),
                        dtype=np.dtype('int32'))
    _dll.openmc_id_map(POINTER(_PlotBase)(plot),
                       img_data.ctypes.data_as(POINTER(c_int32)))
    return img_data
