from ctypes import (c_bool, c_int, c_size_t, c_int32,
                    c_double, Structure, POINTER)

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
    u_span_ : openmc.lib.plot._Position
        Full-width span vector defining the plot's horizontal axis.
    v_span_ : openmc.lib.plot._Position
        Full-height span vector defining the plot's vertical axis.
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
                ('u_span_', _Position),
                ('v_span_', _Position),
                ('width_', _Position),
                ('basis_', c_int),
                ('pixels_', 3*c_size_t),
                ('color_overlaps_', c_bool),
                ('level_', c_int)]

    def __init__(self):
        self.level_ = -1
        self.basis_ = 1
        self.color_overlaps_ = False
        self._update_spans()

    def _update_spans(self):
        if self.basis_ == 1:
            self.u_span_.x = self.width_.x
            self.u_span_.y = 0.0
            self.u_span_.z = 0.0
            self.v_span_.x = 0.0
            self.v_span_.y = self.width_.y
            self.v_span_.z = 0.0
        elif self.basis_ == 2:
            self.u_span_.x = self.width_.x
            self.u_span_.y = 0.0
            self.u_span_.z = 0.0
            self.v_span_.x = 0.0
            self.v_span_.y = 0.0
            self.v_span_.z = self.width_.y
        elif self.basis_ == 3:
            self.u_span_.x = 0.0
            self.u_span_.y = self.width_.x
            self.u_span_.z = 0.0
            self.v_span_.x = 0.0
            self.v_span_.y = 0.0
            self.v_span_.z = self.width_.y

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
        self._update_spans()

    @property
    def height(self):
        return self.width_.y

    @height.setter
    def height(self, height):
        self.width_.y = height
        self._update_spans()

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
            self._update_spans()
            return

        if isinstance(basis, int):
            valid_bases = (1, 2, 3)
            if basis not in valid_bases:
                raise ValueError(f"{basis} is not a valid plot basis.")
            self.basis_ = basis
            self._update_spans()
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
    img_data = np.zeros((plot.v_res, plot.h_res, 3), dtype=np.int32)
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


_dll.openmc_slice_plot.argtypes = [
    POINTER(c_double * 3),   # origin
    POINTER(c_double * 3),   # u_span
    POINTER(c_double * 3),   # v_span
    POINTER(c_size_t * 2),   # pixels
    c_bool,                  # color_overlaps
    c_int,                   # level
    c_int32,                 # filter_index
    POINTER(c_int32),        # geom_data
    POINTER(c_double),       # property_data (can be None)
]
_dll.openmc_slice_plot.restype = c_int
_dll.openmc_slice_plot.errcheck = _error_handler


def slice_plot(origin, width=None, basis='xy', u_span=None, v_span=None,
                pixels=None, color_overlaps=False, level=-1, filter=None,
                include_properties=True):
    """Generate a 2D raster of geometry and property data for plotting.

    Parameters
    ----------
    origin : sequence of float
        Center position of the plot [x, y, z]
    width : sequence of float
        Width of the plot [horizontal, vertical]. Mutually exclusive with
        u_span/v_span.
    basis : {'xy', 'xz', 'yz'} or int
        Plot basis. Ignored if u_span/v_span are provided.
    u_span : sequence of float, optional
        Full-width span vector for the horizontal axis (3 values). Mutually
        exclusive with width.
    v_span : sequence of float, optional
        Full-height span vector for the vertical axis (3 values). Mutually
        exclusive with width.
    pixels : sequence of int
        Number of pixels [horizontal, vertical]
    color_overlaps : bool, optional
        Whether to detect overlapping cells
    level : int, optional
        Universe level (-1 for deepest)
    filter : openmc.lib.Filter, optional
        Filter for bin index lookup
    include_properties : bool, optional
        Whether to compute temperature/density
    Returns
    -------
    geom_data : numpy.ndarray
        Array of shape (v_res, h_res, 3) or (v_res, h_res, 4) with int32 dtype.
        Contains [cell_id, cell_instance, material_id] when no filter is provided,
        or [cell_id, cell_instance, material_id, filter_bin] when a filter is provided.
    property_data : numpy.ndarray or None
        Array of shape (v_res, h_res, 2) with float64 dtype containing
        [temperature, density], or None if include_properties=False
    """
    if pixels is None:
        raise ValueError("pixels must be specified.")
    if len(pixels) != 2:
        raise ValueError("pixels must be a length-2 sequence.")

    if width is not None and (u_span is not None or v_span is not None):
        raise ValueError("width is mutually exclusive with u_span/v_span.")

    if u_span is not None or v_span is not None:
        if u_span is None or v_span is None:
            raise ValueError("Both u_span and v_span must be provided.")
        u_span = np.asarray(u_span, dtype=float)
        v_span = np.asarray(v_span, dtype=float)
        if u_span.shape != (3,) or v_span.shape != (3,):
            raise ValueError("u_span and v_span must be length-3 sequences.")
        u_norm = np.linalg.norm(u_span)
        v_norm = np.linalg.norm(v_span)
        if u_norm == 0.0 or v_norm == 0.0:
            raise ValueError("u_span and v_span must be non-zero vectors.")
        dot = float(np.dot(u_span, v_span))
        ortho_tol = 1.0e-10 * u_norm * v_norm
        if abs(dot) > ortho_tol:
            raise ValueError("u_span and v_span must be orthogonal.")
    else:
        if width is None:
            raise ValueError("width must be provided when u_span/v_span are not set.")
        if len(width) != 2:
            raise ValueError("width must be a length-2 sequence.")
        basis_map = {'xy': 1, 'xz': 2, 'yz': 3}
        if isinstance(basis, str):
            basis = basis.lower()
            if basis not in basis_map:
                raise ValueError(f"{basis} is not a valid plot basis.")
            basis = basis_map[basis]
        elif isinstance(basis, int):
            if basis not in basis_map.values():
                raise ValueError(f"{basis} is not a valid plot basis.")
        else:
            raise ValueError(f"{basis} is not a valid plot basis.")

        if basis == 1:
            u_span = np.array([width[0], 0.0, 0.0], dtype=float)
            v_span = np.array([0.0, width[1], 0.0], dtype=float)
        elif basis == 2:
            u_span = np.array([width[0], 0.0, 0.0], dtype=float)
            v_span = np.array([0.0, 0.0, width[1]], dtype=float)
        else:
            u_span = np.array([0.0, width[0], 0.0], dtype=float)
            v_span = np.array([0.0, 0.0, width[1]], dtype=float)

    origin = np.asarray(origin, dtype=float)
    if origin.shape != (3,):
        raise ValueError("origin must be a length-3 sequence.")

    # Prepare ctypes arrays
    origin_arr = (c_double * 3)(*origin)
    u_span_arr = (c_double * 3)(*u_span)
    v_span_arr = (c_double * 3)(*v_span)
    pixels_arr = (c_size_t * 2)(*pixels)

    # Get internal filter index from filter ID if filter is provided
    if filter is not None:
        filter_index = c_int32()
        _dll.openmc_get_filter_index(filter.id, filter_index)
        filter_index = filter_index.value
    else:
        filter_index = -1

    # Allocate output arrays with dynamic size based on filter
    n_geom_fields = 4 if filter is not None else 3
    geom_data = np.zeros((pixels[1], pixels[0], n_geom_fields), dtype=np.int32)
    if include_properties:
        property_data = np.zeros((pixels[1], pixels[0], 2), dtype=np.float64)
        prop_ptr = property_data.ctypes.data_as(POINTER(c_double))
    else:
        property_data = None
        prop_ptr = None

    _dll.openmc_slice_plot(
        origin_arr,
        u_span_arr,
        v_span_arr,
        pixels_arr,
        color_overlaps,
        level,
        filter_index,
        geom_data.ctypes.data_as(POINTER(c_int32)),
        prop_ptr
    )

    return geom_data, property_data
