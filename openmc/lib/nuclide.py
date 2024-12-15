from collections.abc import Mapping
from ctypes import c_int, c_double, c_char_p, POINTER, c_size_t
from weakref import WeakValueDictionary

from numpy.ctypeslib import ndpointer
import numpy as np

from ..exceptions import DataError, AllocationError
from . import _dll
from .core import _FortranObject
from .error import _error_handler


__all__ = ["Nuclide", "nuclides", "load_nuclide"]

_array_1d_dble = ndpointer(dtype=np.double, ndim=1, flags="CONTIGUOUS")

# Nuclide functions
_dll.openmc_get_nuclide_index.argtypes = [c_char_p, POINTER(c_int)]
_dll.openmc_get_nuclide_index.restype = c_int
_dll.openmc_get_nuclide_index.errcheck = _error_handler
_dll.openmc_load_nuclide.argtypes = [c_char_p, POINTER(c_double), c_int]
_dll.openmc_load_nuclide.restype = c_int
_dll.openmc_load_nuclide.errcheck = _error_handler
_dll.openmc_nuclide_name.argtypes = [c_int, POINTER(c_char_p)]
_dll.openmc_nuclide_name.restype = c_int
_dll.openmc_nuclide_name.errcheck = _error_handler
_dll.openmc_nuclide_collapse_rate.argtypes = [
    c_int,
    c_int,
    c_double,
    _array_1d_dble,
    _array_1d_dble,
    c_int,
    POINTER(c_double),
]
_dll.openmc_nuclide_collapse_rate.restype = c_int
_dll.openmc_nuclide_collapse_rate.errcheck = _error_handler
_dll.nuclides_size.restype = c_size_t


def load_nuclide(name):
    """Load cross section data for a nuclide.

    Parameters
    ----------
    name : str
        Name of the nuclide, e.g. 'U235'

    """
    _dll.openmc_load_nuclide(name.encode(), None, 0)


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

    def collapse_rate(self, MT, temperature, energy, flux):
        """Calculate reaction rate based on group-wise flux distribution

        Parameters
        ----------
        MT : int
            ENDF MT value of the desired reaction
        temperature : float
            Temperature in [K] at which to evaluate cross sections
        energy : iterable of float
            Energy group boundaries in [eV]
        flux : iterable of float
            Flux in each energy group (not normalized per eV)

        Returns
        -------
        float
            Reaction rate

        """
        energy = np.asarray(energy, dtype=float)
        flux = np.asarray(flux, dtype=float)
        xs = c_double()
        _dll.openmc_nuclide_collapse_rate(
            self._index, MT, temperature, energy, flux, len(flux), xs
        )
        return xs.value


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
