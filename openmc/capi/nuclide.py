from collections import Mapping
from ctypes import c_int, c_char_p, POINTER
from weakref import WeakValueDictionary

import numpy as np
from numpy.ctypeslib import as_array

from openmc.capi import _dll, _error_handler


__all__ = ['NuclideView', 'nuclides', 'load_nuclide']

# Nuclide functions
_dll.openmc_get_nuclide.argtypes = [c_char_p, POINTER(c_int)]
_dll.openmc_get_nuclide.restype = c_int
_dll.openmc_get_nuclide.errcheck = _error_handler
_dll.openmc_load_nuclide.argtypes = [c_char_p]
_dll.openmc_load_nuclide.restype = c_int
_dll.openmc_load_nuclide.errcheck = _error_handler
_dll.openmc_nuclide_name.argtypes = [c_int, POINTER(c_char_p)]
_dll.openmc_nuclide_name.restype = c_int
_dll.openmc_nuclide_name.errcheck = _error_handler


def load_nuclide(name):
    """Load cross section data for a nuclide.

    Parameters
    ----------
    name : str
        Name of nuclide, e.g. 'U235'

    """
    _dll.openmc_load_nuclide(name.encode())


class NuclideView(object):
    __instances = WeakValueDictionary()
    def __new__(cls, *args):
        if args not in cls.__instances:
            instance = super().__new__(cls)
            cls.__instances[args] = instance
        return cls.__instances[args]

    def __init__(self, index):
        self._index = index

    @property
    def name(self):
        """Name of nuclide with given index

        Parameter
        ---------
        index : int
            Index in internal nuclides array

        Returns
        -------
        str
            Name of nuclide

        """
        name = c_char_p()
        _dll.openmc_nuclide_name(self._index, name)

        # Find blank in name
        i = 0
        while name.value[i:i+1] != b' ':
            i += 1
        return name.value[:i].decode()


class _NuclideMapping(Mapping):
    def __getitem__(self, key):
        index = c_int()
        _dll.openmc_get_nuclide(key.encode(), index)
        return NuclideView(index)

    def __iter__(self):
        for i in range(len(self)):
            yield NuclideView(i + 1).name

    def __len__(self):
        return c_int.in_dll(_dll, 'n_nuclides').value

nuclides = _NuclideMapping()
