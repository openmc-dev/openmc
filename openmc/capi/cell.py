from collections import Mapping, Iterable
from ctypes import c_int, c_int32, c_double, c_char_p, POINTER
from weakref import WeakValueDictionary

import numpy as np
from numpy.ctypeslib import as_array

from . import _dll
from .core import _View
from .error import _error_handler
from .material import MaterialView

__all__ = ['CellView', 'cells']

# Cell functions
_dll.openmc_cell_get_id.argtypes = [c_int32, POINTER(c_int32)]
_dll.openmc_cell_get_id.restype = c_int
_dll.openmc_cell_get_id.errcheck = _error_handler
_dll.openmc_cell_get_fill.argtypes = [
    c_int32, POINTER(c_int), POINTER(POINTER(c_int32)), POINTER(c_int32)]
_dll.openmc_cell_set_fill.argtypes = [
    c_int32, c_int, c_int32, POINTER(c_int32)]
_dll.openmc_cell_set_fill.restype = c_int
_dll.openmc_cell_set_fill.errcheck = _error_handler
_dll.openmc_cell_set_temperature.argtypes = [
    c_int32, c_double, POINTER(c_int32)]
_dll.openmc_cell_set_temperature.restype = c_int
_dll.openmc_cell_set_temperature.errcheck = _error_handler
_dll.openmc_get_cell_index.argtypes = [c_int32, POINTER(c_int32)]
_dll.openmc_get_cell_index.restype = c_int
_dll.openmc_get_cell_index.errcheck = _error_handler


class CellView(_View):
    """View of a cell.

    This class exposes a cell that is stored internally in the OpenMC solver. To
    obtain a view of a cell with a given ID, use the
    :data:`openmc.capi.nuclides` mapping.

    Parameters
    ----------
    index : int
         Index in the `cells` array.

    Attributes
    ----------
    id : int
        ID of the cell

    """
    __instances = WeakValueDictionary()

    def __new__(cls, *args):
        if args not in cls.__instances:
            instance = super().__new__(cls)
            cls.__instances[args] = instance
        return cls.__instances[args]

    @property
    def id(self):
        cell_id = c_int32()
        _dll.openmc_cell_get_id(self._index, cell_id)
        return cell_id.value

    @property
    def fill(self):
        fill_type = c_int()
        indices = POINTER(c_int32)()
        n = c_int32()
        _dll.openmc_cell_get_fill(self._index, fill_type, indices, n)

        if fill_type.value == 1:
            if n.value > 1:
                return [MaterialView(i) for i in indices[:n.value]]
            else:
                return MaterialView(indices[0])
        else:
            raise NotImplementedError

    @fill.setter
    def fill(self, fill):
        if isinstance(fill, Iterable):
            n = len(fill)
            indices = (c_int*n)(*(m._index for m in fill))
            _dll.openmc_cell_set_fill(self._index, 1, 1, indices)
        elif isinstance(fill, MaterialView):
            materials = [fill]
            indices = (c_int*1)(fill._index)
            _dll.openmc_cell_set_fill(self._index, 1, 1, indices)
        else:
            raise NotImplementedError

    def set_temperature(self, T, instance=None):
        """Set the temperature of a cell

        Parameters
        ----------
        T : float
            Temperature in K
        instance : int or None
            Which instance of the cell

        """
        _dll.openmc_cell_set_temperature(self._index, T, instance)


class _CellMapping(Mapping):
    def __getitem__(self, key):
        index = c_int32()
        _dll.openmc_get_cell_index(key, index)
        return CellView(index.value)

    def __iter__(self):
        for i in range(len(self)):
            yield CellView(i + 1).id

    def __len__(self):
        return c_int32.in_dll(_dll, 'n_cells').value

    def __repr__(self):
        return repr(dict(self))

cells = _CellMapping()
