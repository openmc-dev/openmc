from collections.abc import Mapping
from ctypes import c_int, c_char_p, POINTER, c_size_t
from weakref import WeakValueDictionary

import numpy as np
from numpy.ctypeslib import as_array

from openmc.exceptions import DataError, AllocationError
from . import _dll
from .core import _FortranObject
from .error import _error_handler


__all__ = ['Nuclide', 'nuclides', 'load_nuclide']

# Nuclide functions
_dll.openmc_get_nuclide_index.argtypes = [c_char_p, POINTER(c_int)]
_dll.openmc_get_nuclide_index.restype = c_int
_dll.openmc_get_nuclide_index.errcheck = _error_handler
_dll.openmc_load_nuclide.argtypes = [c_char_p]
_dll.openmc_load_nuclide.restype = c_int
_dll.openmc_load_nuclide.errcheck = _error_handler
_dll.openmc_nuclide_name.argtypes = [c_int, POINTER(c_char_p)]
_dll.openmc_nuclide_name.restype = c_int
_dll.openmc_nuclide_name.errcheck = _error_handler
_dll.nuclides_size.restype = c_size_t


def load_nuclide(name):
    """Load cross section data for a nuclide.

    Parameters
    ----------
    name : str
        Name of the nuclide, e.g. 'U235'

    """
    _dll.openmc_load_nuclide(name.encode())


class Nuclide(_FortranObject):
    """Nuclide stored internally.

    This class exposes a nuclide that is stored internally in the OpenMC
    solver. To obtain a view of a nuclide with a given name, use the
    :data:`openmc.lib.nuclides` mapping.

    Parameters
    ----------
    index : int
         Index in the `nuclides` array.

    Attributes
    ----------
    name : str
        Name of the nuclide, e.g. 'U235'

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
    def name(self):
        name = c_char_p()
        _dll.openmc_nuclide_name(self._index, name)
        return name.value.decode()


class _NuclideMapping(Mapping):
    """Provide mapping from nuclide name to index in nuclides array."""
    def __getitem__(self, key):
        index = c_int()
        try:
            _dll.openmc_get_nuclide_index(key.encode(), index)
        except (DataError, AllocationError) as e:
            # __contains__ expects a KeyError to work correctly
            raise KeyError(str(e))
        return Nuclide(index.value)

    def __iter__(self):
        for i in range(len(self)):
            yield Nuclide(i).name

    def __len__(self):
        return _dll.nuclides_size()

    def __repr__(self):
        return repr(dict(self))

nuclides = _NuclideMapping()
