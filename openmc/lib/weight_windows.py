from collections.abc import Mapping
from ctypes import c_int, c_int32, c_char_p, c_size_t, POINTER
from weakref import WeakValueDictionary

from openmc.exceptions import AllocationError, InvalidIDError
from . import _dll
from .core import _FortranObjectWithID
from .error import _error_handler
import openmc.lib
import openmc

_dll.openmc_extend_weight_windows.argtypes = [c_int32, POINTER(c_int32), POINTER(c_int32)]
_dll.openmc_update_weight_windows.argtypes = 2*[c_int32] + 3*[c_char_p]
_dll.openmc_update_weight_windows.restype = c_int
_dll.openmc_update_weight_windows.errcheck = _error_handler
_dll.openmc_weight_windows_size.restype = c_size_t


class WeightWindows(_FortranObjectWithID):
    """WeightWindows stored internally.

    Parameters
    ----------
    id : int or None
        Unique ID of the weight windows
    new : bool
        When `index` is None, this argument controls whether a new object is
        created or a view of an existing object is returned.
    index : int or None
        Index in the `weight_windows` array.

    """

    __instances = WeakValueDictionary()

    def __new__(cls, id=None, new=True, index=None):
        mapping = weight_windows

        if index is None:
            if new:
                # Determine ID to assign
                if id is None:
                    id = max(mapping, default=0) + 1
                else:
                    if id in mapping:
                        raise AllocationError(f'A weight windows object with ID={id} '
                                              'has already been allocated.')

                index = c_int32()
                _dll.openmc_extend_weight_windows(1, index, None)
                index = index.value
            else:
                index = mapping[id]._index

        if index not in cls.__instances:
            instance = super().__new__(cls)
            instance._index = index
            if id is not None:
                instance.id = id
            cls.__instances[index] = instance

        return cls.__instances[index]

    def update_weight_windows(self, tally, score, value, method):
        """Update weight window values using tally information

        Parameters
        ----------
        tally : openmc.lib.Tally object
            The tally used to update weight window information
        score : str
            Name of the score data used (default is "flux")
        value : str
            Value type used to generate weight windows. One of {'mean', 'rel_err', 'std_dev}.
            (default is 'mean')
        method : str
            Method used for weight window generation. One of {'magic', 'pseudo-source'}

        """

        _dll.openmc_update_weight_windows(tally_idx,
                                          ww_idx,
                                          c_char_p(score.encode()),
                                          c_char_p(value.encode()),
                                          c_char_p(method.encode()))


class _WeightWindowsMapping(Mapping):
    def __getitem__(self, key):
        index = c_int32()
        try:
            _dll.openmc_get_weight_window_index(key, index)
        except(AllocationError, InvalidIDError):
            raise KeyError(str(e))
        return WeightWindows(index=index.value)

    def __iter__(self):
        for i in range(len(self)):
            yield WeightWindows(index=i).id

    def __len__(self):
        return _dll.openmc_weight_windows_size()

    def __repr__(self):
        return repr(dict(self))

    def __delitem__(self):
        raise NotImplementedError("WeightWindows object remove not implemented")

weight_windows = _WeightWindowsMapping()