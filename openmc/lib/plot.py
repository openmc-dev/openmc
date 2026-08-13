from collections.abc import Mapping
from ctypes import (c_bool, c_int, c_size_t, c_int32,
                    c_double, c_uint8, Structure, POINTER)
from weakref import WeakValueDictionary

from ..exceptions import AllocationError, InvalidIDError
from . import _dll
from .core import _FortranObjectWithID
from .error import _error_handler

import numpy as np
import warnings


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


def _extract_slice_data_args(plot):
    """Convert a legacy plot-like object into slice_data keyword arguments."""
    try:
        kwargs = {
            'origin': tuple(plot.origin),
            'width': (plot.width, plot.height),
            'basis': plot.basis,
            'pixels': (plot.h_res, plot.v_res),
            'show_overlaps': getattr(plot, 'color_overlaps', False),
            'level': getattr(plot, 'level', -1),
        }
    except AttributeError as exc:
        raise TypeError(
            "plot must be a legacy plot-like object with origin, width, "
            "height, basis, h_res, and v_res attributes."
        ) from exc
    return kwargs


_dll.openmc_slice_data.argtypes = [
    POINTER(c_double * 3),   # origin
    POINTER(c_double * 3),   # u_span
    POINTER(c_double * 3),   # v_span
    POINTER(c_size_t * 2),   # pixels
    c_bool,                  # show_overlaps
    c_int,                   # level
    c_int32,                 # filter_index
    POINTER(c_int32),        # geom_data
    POINTER(c_double),       # property_data (can be None)
]
_dll.openmc_slice_data.restype = c_int
_dll.openmc_slice_data.errcheck = _error_handler


def slice_data(origin, width=None, basis='xy', u_span=None, v_span=None,
                pixels=None, show_overlaps=False, level=None, filter=None,
                include_properties=True):
    """Generate a 2D raster of geometry and property data for plotting.

    .. versionadded:: 0.16.0

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
    show_overlaps : bool, optional
        Whether to detect overlapping cells
    level : int, optional
        Universe level (None for deepest)
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
    # Set deepest level as default
    if level is None:
        level = -1
    if not isinstance(level, int):
        raise TypeError("level must be an integer.")

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

    _dll.openmc_slice_data(
        origin_arr,
        u_span_arr,
        v_span_arr,
        pixels_arr,
        show_overlaps,
        level,
        filter_index,
        geom_data.ctypes.data_as(POINTER(c_int32)),
        prop_ptr
    )

    return geom_data, property_data


def id_map(plot):
    """Deprecated compatibility wrapper for geometry ID maps.

    This function is kept for compatibility and will be removed in a future
    release. Use `slice_data(..., include_properties=False)` instead.
    """
    warnings.warn(
        "openmc.lib.id_map is deprecated and will be removed in a future "
        "release; use openmc.lib.slice_data(..., include_properties=False).",
        FutureWarning,
    )

    kwargs = _extract_slice_data_args(plot)
    geom_data, _ = slice_data(include_properties=False, **kwargs)
    return geom_data[:, :, :3]


def property_map(plot):
    """Deprecated compatibility wrapper for temperature/density maps.

    This function is kept for compatibility and will be removed in a future
    release. Use `slice_data(..., include_properties=True)` instead.
    """
    warnings.warn(
        "openmc.lib.property_map is deprecated and will be removed in a "
        "future release; use openmc.lib.slice_data(..., "
        "include_properties=True).",
        FutureWarning,
    )

    kwargs = _extract_slice_data_args(plot)
    _, prop_data = slice_data(include_properties=True, **kwargs)
    return prop_data


_dll.openmc_slice_data_overlap_count.argtypes = [POINTER(c_size_t)]
_dll.openmc_slice_data_overlap_count.restype = c_int
_dll.openmc_slice_data_overlap_count.errcheck = _error_handler

_dll.openmc_slice_data_overlap_info.argtypes = [c_size_t, POINTER(c_int32)]
_dll.openmc_slice_data_overlap_info.restype = c_int
_dll.openmc_slice_data_overlap_info.errcheck = _error_handler


# Python wrappings for overlap functions
def slice_data_overlap_count() -> int:
    """Return the number of unique overlaps from the last slice plot.

    .. versionadded:: 0.16.0

    Returns
    -------
    int
        Number of unique overlapping cell pairs detected by the most recent
        :func:`slice_data` call with overlap checking enabled.
    """
    count = c_size_t()
    _dll.openmc_slice_data_overlap_count(count)
    return count.value


def slice_data_overlap_info() -> np.ndarray:
    """Return identifying information for overlaps from the last slice plot.

    .. versionadded:: 0.16.0

    Returns
    -------
    numpy.ndarray
        Array of shape ``(n, 3)`` with int32 dtype, where ``n`` is the number
        of unique overlaps detected by the most recent :func:`slice_data` call
        with overlap checking enabled. Each row contains ``[universe_id,
        cell1_id, cell2_id]``.
    """
    n = slice_data_overlap_count()
    overlap_info = np.empty((n, 3), dtype=np.int32)

    if n > 0:
        _dll.openmc_slice_data_overlap_info(
            n,
            overlap_info.ctypes.data_as(POINTER(c_int32)),
        )
    return overlap_info


_dll.openmc_get_plot_index.argtypes = [c_int32, POINTER(c_int32)]
_dll.openmc_get_plot_index.restype = c_int
_dll.openmc_get_plot_index.errcheck = _error_handler

_dll.openmc_plot_get_id.argtypes = [c_int32, POINTER(c_int32)]
_dll.openmc_plot_get_id.restype = c_int
_dll.openmc_plot_get_id.errcheck = _error_handler

_dll.openmc_plot_set_id.argtypes = [c_int32, c_int32]
_dll.openmc_plot_set_id.restype = c_int
_dll.openmc_plot_set_id.errcheck = _error_handler

_dll.openmc_plots_size.restype = c_size_t

_dll.openmc_solidraytrace_plot_create.argtypes = [POINTER(c_int32)]
_dll.openmc_solidraytrace_plot_create.restype = c_int
_dll.openmc_solidraytrace_plot_create.errcheck = _error_handler

_dll.openmc_solidraytrace_plot_get_pixels.argtypes = [
    c_int32, POINTER(c_int32), POINTER(c_int32)]
_dll.openmc_solidraytrace_plot_get_pixels.restype = c_int
_dll.openmc_solidraytrace_plot_get_pixels.errcheck = _error_handler

_dll.openmc_solidraytrace_plot_set_pixels.argtypes = [c_int32, c_int32, c_int32]
_dll.openmc_solidraytrace_plot_set_pixels.restype = c_int
_dll.openmc_solidraytrace_plot_set_pixels.errcheck = _error_handler

_dll.openmc_solidraytrace_plot_get_color_by.argtypes = [c_int32, POINTER(c_int32)]
_dll.openmc_solidraytrace_plot_get_color_by.restype = c_int
_dll.openmc_solidraytrace_plot_get_color_by.errcheck = _error_handler

_dll.openmc_solidraytrace_plot_set_color_by.argtypes = [c_int32, c_int32]
_dll.openmc_solidraytrace_plot_set_color_by.restype = c_int
_dll.openmc_solidraytrace_plot_set_color_by.errcheck = _error_handler

_dll.openmc_solidraytrace_plot_set_default_colors.argtypes = [c_int32]
_dll.openmc_solidraytrace_plot_set_default_colors.restype = c_int
_dll.openmc_solidraytrace_plot_set_default_colors.errcheck = _error_handler

_dll.openmc_solidraytrace_plot_set_all_opaque.argtypes = [c_int32]
_dll.openmc_solidraytrace_plot_set_all_opaque.restype = c_int
_dll.openmc_solidraytrace_plot_set_all_opaque.errcheck = _error_handler

_dll.openmc_solidraytrace_plot_set_opaque.argtypes = [c_int32, c_int32, c_bool]
_dll.openmc_solidraytrace_plot_set_opaque.restype = c_int
_dll.openmc_solidraytrace_plot_set_opaque.errcheck = _error_handler

_dll.openmc_solidraytrace_plot_set_color.argtypes = [c_int32, c_int32, c_uint8, c_uint8, c_uint8]
_dll.openmc_solidraytrace_plot_set_color.restype = c_int
_dll.openmc_solidraytrace_plot_set_color.errcheck = _error_handler

_dll.openmc_solidraytrace_plot_get_camera_position.argtypes = [
    c_int32, POINTER(c_double), POINTER(c_double), POINTER(c_double)]
_dll.openmc_solidraytrace_plot_get_camera_position.restype = c_int
_dll.openmc_solidraytrace_plot_get_camera_position.errcheck = _error_handler

_dll.openmc_solidraytrace_plot_set_camera_position.argtypes = [c_int32, c_double, c_double, c_double]
_dll.openmc_solidraytrace_plot_set_camera_position.restype = c_int
_dll.openmc_solidraytrace_plot_set_camera_position.errcheck = _error_handler

_dll.openmc_solidraytrace_plot_get_look_at.argtypes = [
    c_int32, POINTER(c_double), POINTER(c_double), POINTER(c_double)]
_dll.openmc_solidraytrace_plot_get_look_at.restype = c_int
_dll.openmc_solidraytrace_plot_get_look_at.errcheck = _error_handler

_dll.openmc_solidraytrace_plot_set_look_at.argtypes = [c_int32, c_double, c_double, c_double]
_dll.openmc_solidraytrace_plot_set_look_at.restype = c_int
_dll.openmc_solidraytrace_plot_set_look_at.errcheck = _error_handler

_dll.openmc_solidraytrace_plot_get_up.argtypes = [
    c_int32, POINTER(c_double), POINTER(c_double), POINTER(c_double)]
_dll.openmc_solidraytrace_plot_get_up.restype = c_int
_dll.openmc_solidraytrace_plot_get_up.errcheck = _error_handler

_dll.openmc_solidraytrace_plot_set_up.argtypes = [c_int32, c_double, c_double, c_double]
_dll.openmc_solidraytrace_plot_set_up.restype = c_int
_dll.openmc_solidraytrace_plot_set_up.errcheck = _error_handler

_dll.openmc_solidraytrace_plot_get_light_position.argtypes = [
    c_int32, POINTER(c_double), POINTER(c_double), POINTER(c_double)]
_dll.openmc_solidraytrace_plot_get_light_position.restype = c_int
_dll.openmc_solidraytrace_plot_get_light_position.errcheck = _error_handler

_dll.openmc_solidraytrace_plot_set_light_position.argtypes = [c_int32, c_double, c_double, c_double]
_dll.openmc_solidraytrace_plot_set_light_position.restype = c_int
_dll.openmc_solidraytrace_plot_set_light_position.errcheck = _error_handler

_dll.openmc_solidraytrace_plot_get_fov.argtypes = [c_int32, POINTER(c_double)]
_dll.openmc_solidraytrace_plot_get_fov.restype = c_int
_dll.openmc_solidraytrace_plot_get_fov.errcheck = _error_handler

_dll.openmc_solidraytrace_plot_set_fov.argtypes = [c_int32, c_double]
_dll.openmc_solidraytrace_plot_set_fov.restype = c_int
_dll.openmc_solidraytrace_plot_set_fov.errcheck = _error_handler

_dll.openmc_solidraytrace_plot_update_view.argtypes = [c_int32]
_dll.openmc_solidraytrace_plot_update_view.restype = c_int
_dll.openmc_solidraytrace_plot_update_view.errcheck = _error_handler

_dll.openmc_solidraytrace_plot_create_image.argtypes = [c_int32, POINTER(c_uint8), c_int32, c_int32]
_dll.openmc_solidraytrace_plot_create_image.restype = c_int
_dll.openmc_solidraytrace_plot_create_image.errcheck = _error_handler

_dll.openmc_solidraytrace_plot_get_color.argtypes = [c_int32, c_int32,
                                             POINTER(c_uint8), POINTER(c_uint8), POINTER(c_uint8)]
_dll.openmc_solidraytrace_plot_get_color.restype = c_int
_dll.openmc_solidraytrace_plot_get_color.errcheck = _error_handler

_dll.openmc_solidraytrace_plot_get_diffuse_fraction.argtypes = [
    c_int32, POINTER(c_double)]
_dll.openmc_solidraytrace_plot_get_diffuse_fraction.restype = c_int
_dll.openmc_solidraytrace_plot_get_diffuse_fraction.errcheck = _error_handler

_dll.openmc_solidraytrace_plot_set_diffuse_fraction.argtypes = [c_int32, c_double]
_dll.openmc_solidraytrace_plot_set_diffuse_fraction.restype = c_int
_dll.openmc_solidraytrace_plot_set_diffuse_fraction.errcheck = _error_handler


class SolidRayTracePlot(_FortranObjectWithID):
    """Solid ray-traced plot stored internally.

    This class exposes a solid ray-traced plot that is stored internally in
    the OpenMC library. To obtain a view of an existing plot with a given ID,
    use the :data:`openmc.lib.plots` mapping.

    .. versionadded:: 0.16.0

    Parameters
    ----------
    uid : int or None
        Unique ID of the plot
    new : bool
        When `index` is None, this argument controls whether a new object is
        created or a view of an existing object is returned.
    index : int or None
        Index in the internal plots array.

    Attributes
    ----------
    id : int
        Unique ID of the plot.
    pixels : tuple of int
        Plot image dimensions as ``(width, height)``.
    color_by : int
        Coloring mode. Use :attr:`COLOR_BY_MATERIAL` or
        :attr:`COLOR_BY_CELL`.
    camera_position : tuple of float
        Camera position as ``(x, y, z)``.
    look_at : tuple of float
        Point the camera is aimed at as ``(x, y, z)``.
    up : tuple of float
        Up direction as ``(x, y, z)``.
    light_position : tuple of float
        Position of the light source as ``(x, y, z)``.
    fov : float
        Horizontal field-of-view angle in degrees.
    diffuse_fraction : float
        Fraction of reflected light treated as diffuse (0 to 1).
    """

    COLOR_BY_MATERIAL = 0
    COLOR_BY_CELL = 1
    __instances = WeakValueDictionary()

    def __new__(cls, uid=None, new=True, index=None):
        mapping = plots
        if index is None:
            if new:
                if uid is not None and uid in mapping:
                    raise AllocationError(
                        f'A plot with ID={uid} has already been allocated.'
                    )
                index = c_int32()
                _dll.openmc_solidraytrace_plot_create(index)
                index = index.value
            else:
                index = mapping[uid]._index

        if index not in cls.__instances:
            instance = super().__new__(cls)
            instance._index = index
            if uid is not None:
                instance.id = uid
            cls.__instances[index] = instance

        return cls.__instances[index]

    def __init__(self, uid=None, new=True, index=None):
        super().__init__(uid, new, index)

    @property
    def id(self):
        plot_id = c_int32()
        _dll.openmc_plot_get_id(self._index, plot_id)
        return plot_id.value

    @id.setter
    def id(self, plot_id):
        _dll.openmc_plot_set_id(self._index, plot_id)

    @staticmethod
    def _get_xyz(getter, index):
        x = c_double()
        y = c_double()
        z = c_double()
        getter(index, x, y, z)
        return (x.value, y.value, z.value)

    @staticmethod
    def _set_xyz(setter, index, xyz):
        x, y, z = xyz
        setter(index, float(x), float(y), float(z))

    @property
    def pixels(self):
        width = c_int32()
        height = c_int32()
        _dll.openmc_solidraytrace_plot_get_pixels(self._index, width, height)
        return (width.value, height.value)

    @pixels.setter
    def pixels(self, pixels):
        width, height = pixels
        _dll.openmc_solidraytrace_plot_set_pixels(
            self._index, int(width), int(height))

    @property
    def color_by(self):
        color_by = c_int32()
        _dll.openmc_solidraytrace_plot_get_color_by(self._index, color_by)
        return color_by.value

    @color_by.setter
    def color_by(self, color_by):
        _dll.openmc_solidraytrace_plot_set_color_by(self._index, int(color_by))

    def set_default_colors(self):
        _dll.openmc_solidraytrace_plot_set_default_colors(self._index)

    def set_all_opaque(self):
        _dll.openmc_solidraytrace_plot_set_all_opaque(self._index)

    def set_visibility(self, domain_id, visible):
        _dll.openmc_solidraytrace_plot_set_opaque(
            self._index, int(domain_id), bool(visible)
        )

    def set_color(self, domain_id, color):
        r, g, b = [int(c) for c in color]
        _dll.openmc_solidraytrace_plot_set_color(
            self._index, int(domain_id), r, g, b)

    @property
    def camera_position(self):
        return self._get_xyz(_dll.openmc_solidraytrace_plot_get_camera_position,
                             self._index)

    @camera_position.setter
    def camera_position(self, position):
        self._set_xyz(_dll.openmc_solidraytrace_plot_set_camera_position,
                      self._index, position)

    @property
    def look_at(self):
        return self._get_xyz(_dll.openmc_solidraytrace_plot_get_look_at,
                             self._index)

    @look_at.setter
    def look_at(self, position):
        self._set_xyz(_dll.openmc_solidraytrace_plot_set_look_at,
                      self._index, position)

    @property
    def up(self):
        return self._get_xyz(_dll.openmc_solidraytrace_plot_get_up, self._index)

    @up.setter
    def up(self, direction):
        self._set_xyz(_dll.openmc_solidraytrace_plot_set_up, self._index,
                      direction)

    @property
    def light_position(self):
        return self._get_xyz(_dll.openmc_solidraytrace_plot_get_light_position,
                             self._index)

    @light_position.setter
    def light_position(self, position):
        self._set_xyz(_dll.openmc_solidraytrace_plot_set_light_position,
                      self._index, position)

    @property
    def fov(self):
        fov = c_double()
        _dll.openmc_solidraytrace_plot_get_fov(self._index, fov)
        return fov.value

    @fov.setter
    def fov(self, fov):
        _dll.openmc_solidraytrace_plot_set_fov(self._index, float(fov))

    def update_view(self):
        _dll.openmc_solidraytrace_plot_update_view(self._index)

    def create_image(self):
        width, height = self.pixels
        image = np.zeros((height, width, 3), dtype=np.uint8)
        _dll.openmc_solidraytrace_plot_create_image(
            self._index,
            image.ctypes.data_as(POINTER(c_uint8)),
            width,
            height
        )
        return image

    def get_color(self, domain_id):
        r = c_uint8()
        g = c_uint8()
        b = c_uint8()
        _dll.openmc_solidraytrace_plot_get_color(
            self._index, int(domain_id), r, g, b)
        return int(r.value), int(g.value), int(b.value)

    @property
    def diffuse_fraction(self):
        value = c_double()
        _dll.openmc_solidraytrace_plot_get_diffuse_fraction(self._index, value)
        return value.value

    @diffuse_fraction.setter
    def diffuse_fraction(self, value):
        _dll.openmc_solidraytrace_plot_set_diffuse_fraction(
            self._index, float(value))

    # Backward-compatible setter aliases
    def set_pixels(self, width, height):
        self.pixels = (width, height)

    def set_color_by(self, color_by):
        self.color_by = color_by

    def set_camera_position(self, x, y, z):
        self.camera_position = (x, y, z)

    def set_look_at(self, x, y, z):
        self.look_at = (x, y, z)

    def set_up(self, x, y, z):
        self.up = (x, y, z)

    def set_light_position(self, x, y, z):
        self.light_position = (x, y, z)

    def set_fov(self, fov):
        self.fov = fov

    def set_diffuse_fraction(self, value):
        self.diffuse_fraction = value


class _PlotMapping(Mapping):
    def __getitem__(self, key):
        index = c_int32()
        try:
            _dll.openmc_get_plot_index(key, index)
        except (AllocationError, InvalidIDError) as e:
            raise KeyError(str(e))
        return SolidRayTracePlot(index=index.value)

    def __iter__(self):
        for i in range(len(self)):
            yield SolidRayTracePlot(index=i).id

    def __len__(self):
        return _dll.openmc_plots_size()

    def __repr__(self):
        return repr(dict(self))


plots = _PlotMapping()
