from collections import Mapping
from ctypes import c_int, c_int32, c_double, c_char_p, POINTER
from weakref import WeakValueDictionary

from numpy.ctypeslib import as_array

from openmc.capi import _dll, _error_handler, NuclideView


__all__ = ['TallyView', 'tallies']

# Tally functions
_dll.openmc_get_tally.argtypes = [c_int32, POINTER(c_int32)]
_dll.openmc_get_tally.restype = c_int
_dll.openmc_get_tally.errcheck = _error_handler
_dll.openmc_tally_id.argtypes = [c_int32, POINTER(c_int32)]
_dll.openmc_tally_id.restype = c_int
_dll.openmc_tally_id.errcheck = _error_handler
_dll.openmc_tally_results.argtypes = [
    c_int32, POINTER(POINTER(c_double)), POINTER(c_int*3)]
_dll.openmc_tally_results.restype = c_int
_dll.openmc_tally_results.errcheck = _error_handler


class TallyView(object):
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
        tally_id = c_int32()
        _dll.openmc_tally_id(self._index, tally_id)
        return tally_id.value

    @property
    def results(self):
        """Get tally results array

        Returns
        -------
        numpy.ndarray
            Array that exposes the internal tally results array

        """
        data = POINTER(c_double)()
        shape = (c_int*3)()
        _dll.openmc_tally_results(self._index, data, shape)
        return as_array(data, tuple(shape[::-1]))


class _TallyMapping(Mapping):
    def __getitem__(self, key):
        index = c_int32()
        _dll.openmc_get_tally(key, index)
        return TallyView(index.value)

    def __iter__(self):
        for i in range(len(self)):
            yield TallyView(i + 1).id

    def __len__(self):
        return c_int32.in_dll(_dll, 'n_tallies').value

tallies = _TallyMapping()
