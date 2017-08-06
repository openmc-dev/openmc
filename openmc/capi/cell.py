from collections import Mapping
from ctypes import c_int, c_int32, c_double, c_char_p, POINTER
from weakref import WeakValueDictionary

import numpy as np

from . import _dll
from .error import _error_handler

__all__ = ['CellView', 'cells']

# Cell functions
_dll.openmc_cell_get_id.argtypes = [c_int32, POINTER(c_int32)]
_dll.openmc_cell_get_id.restype = c_int
_dll.openmc_cell_get_id.errcheck = _error_handler
_dll.openmc_cell_set_temperature.argtypes = [
    c_int32, c_double, POINTER(c_int32)]
_dll.openmc_cell_set_temperature.restype = c_int
_dll.openmc_cell_set_temperature.errcheck = _error_handler
_dll.openmc_get_cell.argtypes = [c_int32, POINTER(c_int32)]
_dll.openmc_get_cell.restype = c_int
_dll.openmc_get_cell.errcheck = _error_handler


class CellView(object):
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

    def __init__(self, index):
        self._index = index

    @property
    def id(self):
        cell_id = c_int32()
        _dll.openmc_cell_get_id(self._index, cell_id)
        return cell_id.value

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
        _dll.openmc_get_cell(key, index)
        return CellView(index.value)

    def __iter__(self):
        for i in range(len(self)):
            yield CellView(i + 1).id

    def __len__(self):
        return c_int32.in_dll(_dll, 'n_cells').value

cells = _CellMapping()
