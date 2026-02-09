from ctypes import (c_bool, c_int, c_size_t, c_int32,
                    c_double, c_uint8, c_void_p,
                    Structure, POINTER, byref)

from . import _dll
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

    def __getitem__(self, idx):
        if idx == 0:
            return self.x
        elif idx == 1:
            return self.y
        elif idx == 2:
            return self.z
        else:
            raise IndexError(f"{idx} index is invalid for _Position")

    def __setitem__(self, idx, val):
        if idx == 0:
            self.x = val
        elif idx == 1:
            self.y = val
        elif idx == 2:
            self.z = val
        else:
            raise IndexError(f"{idx} index is invalid for _Position")

    def __repr__(self):
        return f"({self.x}, {self.y}, {self.z})"


class _PlotBase(Structure):
    """A structure defining a 2-D geometry slice with underlying c-types

    C-Type Attributes
    -----------------
    origin_ : openmc.lib.plot._Position
        A position defining the origin of the plot.
    width_ : openmc.lib.plot._Position
        The width of the plot along the x, y, and z axes, respectively
    basis_ : c_int
        The axes basis of the plot view.
    pixels_ : c_size_t[3]
        The resolution of the plot in the horizontal and vertical dimensions
    color_overlaps_ : c_bool
        Whether to assign unique IDs (-3) to overlapping regions.
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
    h_res : int
        The horizontal resolution of the plot in pixels
    v_res : int
        The vertical resolution of the plot in pixels
    level : int
        The universe level for the plot (default: -1 -> all universes shown)
    """
    _fields_ = [('origin_', _Position),
                ('width_', _Position),
                ('basis_', c_int),
                ('pixels_', 3*c_size_t),
                ('color_overlaps_', c_bool),
                ('level_', c_int)]

    def __init__(self):
        self.level_ = -1
        self.basis_ = 1
        self.color_overlaps_ = False

    @property
    def origin(self):
        return self.origin_

    @origin.setter
    def origin(self, origin):
        self.origin_.x = origin[0]
        self.origin_.y = origin[1]
        self.origin_.z = origin[2]

    @property
    def width(self):
        return self.width_.x

    @width.setter
    def width(self, width):
        self.width_.x = width

    @property
    def height(self):
        return self.width_.y

    @height.setter
    def height(self, height):
        self.width_.y = height

    @property
    def basis(self):
        if self.basis_ == 1:
            return 'xy'
        elif self.basis_ == 2:
            return 'xz'
        elif self.basis_ == 3:
            return 'yz'

        raise ValueError(f"Plot basis {self.basis_} is invalid")

    @basis.setter
    def basis(self, basis):
        if isinstance(basis, str):
            valid_bases = ('xy', 'xz', 'yz')
            basis = basis.lower()
            if basis not in valid_bases:
                raise ValueError(f"{basis} is not a valid plot basis.")

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
                raise ValueError(f"{basis} is not a valid plot basis.")
            self.basis_ = basis
            return

        raise ValueError(f"{basis} of type {type(basis)} is an invalid plot basis")

    @property
    def h_res(self):
        return self.pixels_[0]

    @h_res.setter
    def h_res(self, h_res):
        self.pixels_[0] = h_res

    @property
    def v_res(self):
        return self.pixels_[1]

    @v_res.setter
    def v_res(self, v_res):
        self.pixels_[1] = v_res

    @property
    def level(self):
        return int(self.level_)

    @level.setter
    def level(self, level):
        self.level_ = level

    @property
    def color_overlaps(self):
        return self.color_overlaps_

    @color_overlaps.setter
    def color_overlaps(self, color_overlaps):
        self.color_overlaps_ = color_overlaps

    def __repr__(self):
        out_str = ["-----",
                   "Plot:",
                   "-----",
                   f"Origin: {self.origin}",
                   f"Width: {self.width}",
                   f"Height: {self.height}",
                   f"Basis: {self.basis}",
                   f"HRes: {self.h_res}",
                   f"VRes: {self.v_res}",
                   f"Color Overlaps: {self.color_overlaps}",
                   f"Level: {self.level}"]
        return '\n'.join(out_str)


_dll.openmc_id_map.argtypes = [POINTER(_PlotBase), POINTER(c_int32)]
_dll.openmc_id_map.restype = c_int
_dll.openmc_id_map.errcheck = _error_handler


def id_map(plot):
    """
    Generate a 2-D map of cell and material IDs. Used for in-memory image
    generation.

    Parameters
    ----------
    plot : openmc.lib.plot._PlotBase
        Object describing the slice of the model to be generated

    Returns
    -------
    id_map : numpy.ndarray
        A NumPy array with shape (vertical pixels, horizontal pixels, 3) of
        OpenMC property ids with dtype int32. The last dimension of the array
        contains, in order, cell IDs, cell instances, and material IDs.

    """
    img_data = np.zeros((plot.v_res, plot.h_res, 3),
                        dtype=np.dtype('int32'))
    _dll.openmc_id_map(plot, img_data.ctypes.data_as(POINTER(c_int32)))
    return img_data


_dll.openmc_property_map.argtypes = [POINTER(_PlotBase), POINTER(c_double)]
_dll.openmc_property_map.restype = c_int
_dll.openmc_property_map.errcheck = _error_handler


def property_map(plot):
    """
    Generate a 2-D map of cell temperatures and material densities. Used for
    in-memory image generation.

    Parameters
    ----------
    plot : openmc.lib.plot._PlotBase
        Object describing the slice of the model to be generated

    Returns
    -------
    property_map : numpy.ndarray
        A NumPy array with shape (vertical pixels, horizontal pixels, 2) of
        OpenMC property ids with dtype float

    """
    prop_data = np.zeros((plot.v_res, plot.h_res, 2))
    _dll.openmc_property_map(plot, prop_data.ctypes.data_as(POINTER(c_double)))
    return prop_data

_dll.openmc_solidraytrace_plot_create.argtypes = [POINTER(c_void_p)]
_dll.openmc_solidraytrace_plot_create.restype = c_int
_dll.openmc_solidraytrace_plot_create.errcheck = _error_handler

_dll.openmc_solidraytrace_plot_free.argtypes = [c_void_p]
_dll.openmc_solidraytrace_plot_free.restype = c_int
_dll.openmc_solidraytrace_plot_free.errcheck = _error_handler

_dll.openmc_solidraytrace_plot_set_pixels.argtypes = [c_void_p, c_int32, c_int32]
_dll.openmc_solidraytrace_plot_set_pixels.restype = c_int
_dll.openmc_solidraytrace_plot_set_pixels.errcheck = _error_handler

_dll.openmc_solidraytrace_plot_set_color_by.argtypes = [c_void_p, c_int32]
_dll.openmc_solidraytrace_plot_set_color_by.restype = c_int
_dll.openmc_solidraytrace_plot_set_color_by.errcheck = _error_handler

_dll.openmc_solidraytrace_plot_set_default_colors.argtypes = [c_void_p]
_dll.openmc_solidraytrace_plot_set_default_colors.restype = c_int
_dll.openmc_solidraytrace_plot_set_default_colors.errcheck = _error_handler

_dll.openmc_solidraytrace_plot_set_all_opaque.argtypes = [c_void_p]
_dll.openmc_solidraytrace_plot_set_all_opaque.restype = c_int
_dll.openmc_solidraytrace_plot_set_all_opaque.errcheck = _error_handler

_dll.openmc_solidraytrace_plot_set_visibility.argtypes = [c_void_p, c_int32, c_bool]
_dll.openmc_solidraytrace_plot_set_visibility.restype = c_int
_dll.openmc_solidraytrace_plot_set_visibility.errcheck = _error_handler

_dll.openmc_solidraytrace_plot_set_color.argtypes = [c_void_p, c_int32, c_uint8, c_uint8, c_uint8]
_dll.openmc_solidraytrace_plot_set_color.restype = c_int
_dll.openmc_solidraytrace_plot_set_color.errcheck = _error_handler

_dll.openmc_solidraytrace_plot_set_camera_position.argtypes = [c_void_p, c_double, c_double, c_double]
_dll.openmc_solidraytrace_plot_set_camera_position.restype = c_int
_dll.openmc_solidraytrace_plot_set_camera_position.errcheck = _error_handler

_dll.openmc_solidraytrace_plot_set_look_at.argtypes = [c_void_p, c_double, c_double, c_double]
_dll.openmc_solidraytrace_plot_set_look_at.restype = c_int
_dll.openmc_solidraytrace_plot_set_look_at.errcheck = _error_handler

_dll.openmc_solidraytrace_plot_set_up.argtypes = [c_void_p, c_double, c_double, c_double]
_dll.openmc_solidraytrace_plot_set_up.restype = c_int
_dll.openmc_solidraytrace_plot_set_up.errcheck = _error_handler

_dll.openmc_solidraytrace_plot_set_light_position.argtypes = [c_void_p, c_double, c_double, c_double]
_dll.openmc_solidraytrace_plot_set_light_position.restype = c_int
_dll.openmc_solidraytrace_plot_set_light_position.errcheck = _error_handler

_dll.openmc_solidraytrace_plot_set_fov.argtypes = [c_void_p, c_double]
_dll.openmc_solidraytrace_plot_set_fov.restype = c_int
_dll.openmc_solidraytrace_plot_set_fov.errcheck = _error_handler

_dll.openmc_solidraytrace_plot_update_view.argtypes = [c_void_p]
_dll.openmc_solidraytrace_plot_update_view.restype = c_int
_dll.openmc_solidraytrace_plot_update_view.errcheck = _error_handler

_dll.openmc_solidraytrace_plot_create_image.argtypes = [c_void_p, POINTER(c_uint8), c_int32, c_int32]
_dll.openmc_solidraytrace_plot_create_image.restype = c_int
_dll.openmc_solidraytrace_plot_create_image.errcheck = _error_handler

_dll.openmc_solidraytrace_plot_get_color.argtypes = [c_void_p, c_int32,
                                             POINTER(c_uint8), POINTER(c_uint8), POINTER(c_uint8)]
_dll.openmc_solidraytrace_plot_get_color.restype = c_int
_dll.openmc_solidraytrace_plot_get_color.errcheck = _error_handler

_dll.openmc_solidraytrace_plot_set_diffuse_fraction.argtypes = [c_void_p, c_double]
_dll.openmc_solidraytrace_plot_set_diffuse_fraction.restype = c_int
_dll.openmc_solidraytrace_plot_set_diffuse_fraction.errcheck = _error_handler


class SolidRayTracePlot:
    """C API wrapper for ray-traced solid plot generation."""
    COLOR_BY_MATERIAL = 0
    COLOR_BY_CELL = 1

    def __init__(self):
        self._ptr = c_void_p()
        _dll.openmc_solidraytrace_plot_create(byref(self._ptr))
        self._width = None
        self._height = None

    def __del__(self):
        if getattr(self, "_ptr", None):
            try:
                _dll.openmc_solidraytrace_plot_free(self._ptr)
            except Exception:
                # Avoid raising from garbage collection
                pass
            self._ptr = None

    def set_pixels(self, width, height):
        _dll.openmc_solidraytrace_plot_set_pixels(self._ptr, int(width), int(height))
        self._width = int(width)
        self._height = int(height)

    def set_color_by(self, color_by):
        _dll.openmc_solidraytrace_plot_set_color_by(self._ptr, int(color_by))

    def set_default_colors(self):
        _dll.openmc_solidraytrace_plot_set_default_colors(self._ptr)

    def set_all_opaque(self):
        _dll.openmc_solidraytrace_plot_set_all_opaque(self._ptr)

    def set_visibility(self, domain_id, visible):
        _dll.openmc_solidraytrace_plot_set_visibility(self._ptr, int(domain_id), bool(visible))

    def set_color(self, domain_id, color):
        r, g, b = [int(c) for c in color]
        _dll.openmc_solidraytrace_plot_set_color(self._ptr, int(domain_id), r, g, b)

    def set_camera_position(self, x, y, z):
        _dll.openmc_solidraytrace_plot_set_camera_position(self._ptr, float(x), float(y), float(z))

    def set_look_at(self, x, y, z):
        _dll.openmc_solidraytrace_plot_set_look_at(self._ptr, float(x), float(y), float(z))

    def set_up(self, x, y, z):
        _dll.openmc_solidraytrace_plot_set_up(self._ptr, float(x), float(y), float(z))

    def set_light_position(self, x, y, z):
        _dll.openmc_solidraytrace_plot_set_light_position(self._ptr, float(x), float(y), float(z))

    def set_fov(self, fov):
        _dll.openmc_solidraytrace_plot_set_fov(self._ptr, float(fov))

    def update_view(self):
        _dll.openmc_solidraytrace_plot_update_view(self._ptr)

    def create_image(self):
        if self._width is None or self._height is None:
            raise RuntimeError("SolidRayTracePlot pixels must be set before create_image")
        image = np.zeros((self._height, self._width, 3), dtype=np.uint8)
        _dll.openmc_solidraytrace_plot_create_image(
            self._ptr,
            image.ctypes.data_as(POINTER(c_uint8)),
            self._width,
            self._height
        )
        return image

    def get_color(self, domain_id):
        r = c_uint8()
        g = c_uint8()
        b = c_uint8()
        _dll.openmc_solidraytrace_plot_get_color(self._ptr, int(domain_id), r, g, b)
        return int(r.value), int(g.value), int(b.value)

    def set_diffuse_fraction(self, value):
        _dll.openmc_solidraytrace_plot_set_diffuse_fraction(self._ptr, float(value))
