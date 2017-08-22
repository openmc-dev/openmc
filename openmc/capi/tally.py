from collections import Mapping
from ctypes import c_int, c_int32, c_double, c_char_p, POINTER
from weakref import WeakValueDictionary

from numpy.ctypeslib import as_array

from . import _dll, NuclideView
from .error import _error_handler
from .filter import _get_filter


__all__ = ['TallyView', 'tallies']

# Tally functions
_dll.openmc_get_tally.argtypes = [c_int32, POINTER(c_int32)]
_dll.openmc_get_tally.restype = c_int
_dll.openmc_get_tally.errcheck = _error_handler
_dll.openmc_extend_tallies.argtypes = [c_int32, POINTER(c_int32), POINTER(c_int32)]
_dll.openmc_extend_tallies.restype = c_int
_dll.openmc_extend_tallies.errcheck = _error_handler
_dll.openmc_tally_get_id.argtypes = [c_int32, POINTER(c_int32)]
_dll.openmc_tally_get_id.restype = c_int
_dll.openmc_tally_get_id.errcheck = _error_handler
_dll.openmc_tally_get_filters.argtypes = [
    c_int32, POINTER(POINTER(c_int32)), POINTER(c_int)]
_dll.openmc_tally_get_filters.restype = c_int
_dll.openmc_tally_get_filters.errcheck = _error_handler
_dll.openmc_tally_get_nuclides.argtypes = [
    c_int32, POINTER(POINTER(c_int)), POINTER(c_int)]
_dll.openmc_tally_get_nuclides.restype = c_int
_dll.openmc_tally_get_nuclides.errcheck = _error_handler
_dll.openmc_tally_results.argtypes = [
    c_int32, POINTER(POINTER(c_double)), POINTER(c_int*3)]
_dll.openmc_tally_results.restype = c_int
_dll.openmc_tally_results.errcheck = _error_handler
_dll.openmc_tally_set_filters.argtypes = [c_int32, c_int, POINTER(c_int32)]
_dll.openmc_tally_set_filters.restype = c_int
_dll.openmc_tally_set_filters.errcheck = _error_handler
_dll.openmc_tally_set_id.argtypes = [c_int32, c_int32]
_dll.openmc_tally_set_id.restype = c_int
_dll.openmc_tally_set_id.errcheck = _error_handler
_dll.openmc_tally_set_nuclides.argtypes = [c_int32, c_int, POINTER(c_char_p)]
_dll.openmc_tally_set_nuclides.restype = c_int
_dll.openmc_tally_set_nuclides.errcheck = _error_handler
_dll.openmc_tally_set_scores.argtypes = [c_int32, c_int, POINTER(c_char_p)]
_dll.openmc_tally_set_scores.restype = c_int
_dll.openmc_tally_set_scores.errcheck = _error_handler
_dll.openmc_tally_set_type.argtypes = [c_int32, c_char_p]
_dll.openmc_tally_set_type.restype = c_int
_dll.openmc_tally_set_type.errcheck = _error_handler


class TallyView(object):
    """View of a tally.

    This class exposes a tally that is stored internally in the OpenMC
    solver. To obtain a view of a tally with a given ID, use the
    :data:`openmc.capi.tallies` mapping.

    Parameters
    ----------
    index : int
         Index in the `tallies` array.

    Attributes
    ----------
    id : int
        ID of the tally
    filters : list
        List of views to tally filters
    nuclides : list of str
        List of nuclides to score results for
    results : numpy.ndarray
        Array of tally results

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
        tally_id = c_int32()
        _dll.openmc_tally_get_id(self._index, tally_id)
        return tally_id.value

    @id.setter
    def id(self, tally_id):
        _dll.openmc_tally_set_id(self._index, tally_id)

    @property
    def filters(self):
        filt_idx = POINTER(c_int32)()
        n = c_int()
        _dll.openmc_tally_get_filters(self._index, filt_idx, n)
        return [_get_filter(filt_idx[i]) for i in range(n.value)]

    @property
    def nuclides(self):
        nucs = POINTER(c_int)()
        n = c_int()
        _dll.openmc_tally_get_nuclides(self._index, nucs, n)
        return [NuclideView(nucs[i]).name if nucs[i] > 0 else 'total'
                for i in range(n.value)]

    @property
    def results(self):
        data = POINTER(c_double)()
        shape = (c_int*3)()
        _dll.openmc_tally_results(self._index, data, shape)
        return as_array(data, tuple(shape[::-1]))

    @filters.setter
    def filters(self, filters):
        # Get filter indices as int32_t[]
        n = len(filters)
        indices = (c_int32*n)(*(f._index for f in filters))

        _dll.openmc_tally_set_filters(self._index, n, indices)

    @nuclides.setter
    def nuclides(self, nuclides):
        nucs = (c_char_p * len(nuclides))()
        nucs[:] = [x.encode() for x in nuclides]
        _dll.openmc_tally_set_nuclides(self._index, len(nuclides), nucs)

    @property
    def scores(self):
        pass

    @scores.setter
    def scores(self, scores):
        scores_ = (c_char_p * len(scores))()
        scores_[:] = [x.encode() for x in scores]
        _dll.openmc_tally_set_scores(self._index, len(scores), scores_)

    @classmethod
    def new(cls):
        index = c_int32()
        _dll.openmc_extend_tallies(1, index, None)
        _dll.openmc_tally_set_type(index, b'generic')
        return cls(index.value)


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
