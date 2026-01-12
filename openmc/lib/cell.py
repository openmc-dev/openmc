import sys

from collections.abc import Mapping, Iterable
from ctypes import c_int, c_int32, c_double, c_char_p, POINTER, c_bool, c_size_t
from weakref import WeakValueDictionary

import numpy as np

from ..exceptions import AllocationError, InvalidIDError
from . import _dll
from .core import _FortranObjectWithID
from .error import _error_handler
from .material import Material
from ..bounding_box import BoundingBox


__all__ = ['Cell', 'cells']

# Cell functions
_dll.openmc_extend_cells.argtypes = [c_int32, POINTER(c_int32), POINTER(c_int32)]
_dll.openmc_extend_cells.restype = c_int
_dll.openmc_extend_cells.errcheck = _error_handler
_dll.openmc_cell_get_id.argtypes = [c_int32, POINTER(c_int32)]
_dll.openmc_cell_get_id.restype = c_int
_dll.openmc_cell_get_id.errcheck = _error_handler
_dll.openmc_cell_get_fill.argtypes = [
    c_int32, POINTER(c_int), POINTER(POINTER(c_int32)), POINTER(c_int32)]
_dll.openmc_cell_get_fill.restype = c_int
_dll.openmc_cell_get_fill.errcheck = _error_handler
_dll.openmc_cell_get_num_instances.argtypes = [c_int32, POINTER(c_int32)]
_dll.openmc_cell_get_num_instances.restype = c_int
_dll.openmc_cell_get_num_instances.errcheck = _error_handler
_dll.openmc_cell_get_temperature.argtypes = [
    c_int32, POINTER(c_int32), POINTER(c_double)]
_dll.openmc_cell_get_temperature.restype = c_int
_dll.openmc_cell_get_temperature.errcheck = _error_handler
_dll.openmc_cell_get_density.argtypes = [
    c_int32, POINTER(c_int32), POINTER(c_double)]
_dll.openmc_cell_get_density.restype = c_int
_dll.openmc_cell_get_density.errcheck = _error_handler
_dll.openmc_cell_get_name.argtypes = [c_int32, POINTER(c_char_p)]
_dll.openmc_cell_get_name.restype = c_int
_dll.openmc_cell_get_name.errcheck = _error_handler
_dll.openmc_cell_get_translation.argtypes = [c_int32, POINTER(c_double)]
_dll.openmc_cell_get_translation.restype = c_int
_dll.openmc_cell_get_translation.errcheck = _error_handler
_dll.openmc_cell_get_rotation.argtypes = [c_int32, POINTER(c_double),
    POINTER(c_size_t)]
_dll.openmc_cell_get_rotation.restype = c_int
_dll.openmc_cell_get_rotation.errcheck = _error_handler
_dll.openmc_cell_set_name.argtypes = [c_int32, c_char_p]
_dll.openmc_cell_set_name.restype = c_int
_dll.openmc_cell_set_name.errcheck = _error_handler
_dll.openmc_cell_set_fill.argtypes = [
    c_int32, c_int, c_int32, POINTER(c_int32)]
_dll.openmc_cell_set_fill.restype = c_int
_dll.openmc_cell_set_fill.errcheck = _error_handler
_dll.openmc_cell_set_id.argtypes = [c_int32, c_int32]
_dll.openmc_cell_set_id.restype = c_int
_dll.openmc_cell_set_id.errcheck = _error_handler
_dll.openmc_cell_set_temperature.argtypes = [
    c_int32, c_double, POINTER(c_int32), c_bool]
_dll.openmc_cell_set_temperature.restype = c_int
_dll.openmc_cell_set_temperature.errcheck = _error_handler
_dll.openmc_cell_set_density.argtypes = [
    c_int32, c_double, POINTER(c_int32), c_bool]
_dll.openmc_cell_set_density.restype = c_int
_dll.openmc_cell_set_density.errcheck = _error_handler
_dll.openmc_cell_set_translation.argtypes = [c_int32, POINTER(c_double)]
_dll.openmc_cell_set_translation.restype = c_int
_dll.openmc_cell_set_translation.errcheck = _error_handler
_dll.openmc_cell_set_rotation.argtypes = [
    c_int32, POINTER(c_double), c_size_t]
_dll.openmc_cell_set_rotation.restype = c_int
_dll.openmc_cell_set_rotation.errcheck = _error_handler
_dll.openmc_get_cell_index.argtypes = [c_int32, POINTER(c_int32)]
_dll.openmc_get_cell_index.restype = c_int
_dll.openmc_get_cell_index.errcheck = _error_handler
_dll.cells_size.restype = c_int
_dll.openmc_cell_bounding_box.argtypes = [c_int,
                                          POINTER(c_double),
                                          POINTER(c_double)]
_dll.openmc_cell_bounding_box.restype = c_int
_dll.openmc_cell_bounding_box.errcheck = _error_handler


class Cell(_FortranObjectWithID):
    """Cell stored internally.

    This class exposes a cell that is stored internally in the OpenMC
    library. To obtain a view of a cell with a given ID, use the
    :data:`openmc.lib.cells` mapping.

    Parameters
    ----------
    uid : int or None
         Unique ID of the cell
    new : bool
         When `index` is None, this argument controls whether a new object is
         created or a view to an existing object is returned.
    index : int
         Index in the `cells` array.

    Attributes
    ----------
    id : int
        ID of the cell
    fill : openmc.lib.Material or list of openmc.lib.Material
        Indicates what the region of space is filled with
    name : str
        Name of the cell
    num_instances : int
        Number of unique cell instances
    bounding_box : openmc.BoundingBox
        Axis-aligned bounding box of the cell
    translation : Iterable of float
        3-D coordinates of the translation vector
    rotation : Iterable of float
        The rotation matrix or angles of the universe filling the cell. This
        can either be a fully specified 3 x 3 rotation matrix or an Iterable
        of length 3 with the angles in degrees about the x, y, and z axes,
        respectively.

    """
    __instances = WeakValueDictionary()

    def __new__(cls, uid=None, new=True, index=None):
        mapping = cells
        if index is None:
            if new:
                # Determine ID to assign
                if uid is None:
                    uid = max(mapping, default=0) + 1
                else:
                    if uid in mapping:
                        raise AllocationError('A cell with ID={} has already '
                                              'been allocated.'.format(uid))

                index = c_int32()
                _dll.openmc_extend_cells(1, index, None)
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

    @property
    def id(self):
        cell_id = c_int32()
        _dll.openmc_cell_get_id(self._index, cell_id)
        return cell_id.value

    @id.setter
    def id(self, cell_id):
        _dll.openmc_cell_set_id(self._index, cell_id)

    @property
    def name(self):
        name = c_char_p()
        _dll.openmc_cell_get_name(self._index, name)
        return name.value.decode()

    @name.setter
    def name(self, name):
        name_ptr = c_char_p(name.encode())
        _dll.openmc_cell_set_name(self._index, name_ptr)

    @property
    def fill(self):
        fill_type = c_int()
        indices = POINTER(c_int32)()
        n = c_int32()
        _dll.openmc_cell_get_fill(self._index, fill_type, indices, n)
        if fill_type.value == 0:
            if n.value > 1:
                return [Material(index=i) for i in indices[:n.value]]
            else:
                index = indices[0]
                return Material(index=index)
        else:
            raise NotImplementedError

    @fill.setter
    def fill(self, fill):
        if isinstance(fill, Iterable):
            n = len(fill)
            indices = (c_int32*n)(*(m._index if m is not None else -1
                                    for m in fill))
            _dll.openmc_cell_set_fill(self._index, 0, n, indices)
        elif isinstance(fill, Material):
            indices = (c_int32*1)(fill._index)
            _dll.openmc_cell_set_fill(self._index, 0, 1, indices)
        elif fill is None:
            indices = (c_int32*1)(-1)
            _dll.openmc_cell_set_fill(self._index, 0, 1, indices)

    @property
    def num_instances(self):
        n = c_int32()
        _dll.openmc_cell_get_num_instances(self._index, n)
        return n.value

    def get_temperature(self, instance=None):
        """Get the temperature of a cell

        Parameters
        ----------
        instance: int or None
            Which instance of the cell

        """

        if instance is not None:
            instance = c_int32(instance)

        T = c_double()
        _dll.openmc_cell_get_temperature(self._index, instance, T)
        return T.value

    def set_temperature(self, T, instance=None, set_contained=False):
        """Set the temperature of a cell

        Parameters
        ----------
        T : float
            Temperature in K
        instance : int or None
            Which instance of the cell
        set_contained: bool
            If cell is not filled by a material, whether to set the temperature of
            all filled cells

        """

        if instance is not None:
            instance = c_int32(instance)

        _dll.openmc_cell_set_temperature(self._index, T, instance, set_contained)

    def get_density(self, instance: int | None = None):
        """Get the density of a cell in [g/cm3]

        Parameters
        ----------
        instance : int or None
            Which instance of the cell

        """

        if instance is not None:
            instance = c_int32(instance)

        rho = c_double()
        _dll.openmc_cell_get_density(self._index, instance, rho)
        return rho.value

    def set_density(self, rho: float, instance: int | None = None,
                    set_contained: bool = False):
        """Set the density of a cell

        Parameters
        ----------
        rho : float
            Density of the cell in [g/cm3]
        instance : int or None
            Which instance of the cell
        set_contained : bool
            If cell is not filled by a material, whether to set the density
            of all filled cells

        """

        if instance is not None:
            instance = c_int32(instance)

        _dll.openmc_cell_set_density(self._index, rho, instance, set_contained)

    @property
    def translation(self):
        translation = np.zeros(3)
        _dll.openmc_cell_get_translation(
            self._index, translation.ctypes.data_as(POINTER(c_double)))
        return translation

    @translation.setter
    def translation(self, translation_vec):
        vector = np.asarray(translation_vec, dtype=float)
        _dll.openmc_cell_set_translation(
            self._index, vector.ctypes.data_as(POINTER(c_double)))

    @property
    def rotation(self):
        rotation_data = np.zeros(12)
        rot_size = c_size_t()

        _dll.openmc_cell_get_rotation(
            self._index, rotation_data.ctypes.data_as(POINTER(c_double)),
            rot_size)
        rot_size = rot_size.value

        if rot_size == 9:
            return rotation_data[:rot_size].shape(3, 3)
        elif rot_size in (0, 12):
            # If size is 0, rotation_data[9:] will be zeros. This indicates no
            # rotation and is the most straightforward way to always return
            # an iterable of floats
            return rotation_data[9:]
        else:
            raise ValueError(
                f'Invalid size of rotation matrix: {rot_size}')

    @rotation.setter
    def rotation(self, rotation_data):
        flat_rotation = np.asarray(rotation_data, dtype=float).flatten()

        _dll.openmc_cell_set_rotation(
            self._index, flat_rotation.ctypes.data_as(POINTER(c_double)),
            c_size_t(len(flat_rotation)))

    @property
    def bounding_box(self):
        inf = sys.float_info.max
        llc = np.zeros(3)
        urc = np.zeros(3)
        _dll.openmc_cell_bounding_box(self._index,
                                 llc.ctypes.data_as(POINTER(c_double)),
                                 urc.ctypes.data_as(POINTER(c_double)))
        llc[llc == inf] = np.inf
        urc[urc == inf] = np.inf
        llc[llc == -inf] = -np.inf
        urc[urc == -inf] = -np.inf

        return BoundingBox(llc, urc)


class _CellMapping(Mapping):
    def __getitem__(self, key):
        index = c_int32()
        try:
            _dll.openmc_get_cell_index(key, index)
        except (AllocationError, InvalidIDError) as e:
            # __contains__ expects a KeyError to work correctly
            raise KeyError(str(e))
        return Cell(index=index.value)

    def __iter__(self):
        for i in range(len(self)):
            yield Cell(index=i).id

    def __len__(self):
        return _dll.cells_size()

    def __repr__(self):
        return repr(dict(self))

cells = _CellMapping()
