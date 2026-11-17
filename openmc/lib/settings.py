from ctypes import (
    byref, c_bool, c_char_p, c_double, c_int, c_int32, c_int64, POINTER
)

from . import _dll
from .error import _error_handler

_RUN_MODES = {1: 'fixed source',
              2: 'eigenvalue',
              3: 'plot',
              4: 'particle restart',
              5: 'volume'}

_dll.openmc_set_seed.argtypes = [c_int64]
_dll.openmc_get_seed.restype = c_int64
_dll.openmc_set_stride.argtypes = [c_int64]
_dll.openmc_get_stride.restype = c_int64
_dll.openmc_get_n_batches.argtypes = [POINTER(c_int), c_bool]
_dll.openmc_get_n_batches.restype = c_int
_dll.openmc_get_n_batches.errcheck = _error_handler
_dll.openmc_set_n_batches.argtypes = [c_int32, c_bool, c_bool]
_dll.openmc_set_n_batches.restype = c_int
_dll.openmc_set_n_batches.errcheck = _error_handler


_SETTING_ACCESSORS = {}
for _suffix, _ctype in (
    ('bool', c_bool),
    ('int32', c_int32),
    ('int64', c_int64),
    ('double', c_double),
):
    _getter = getattr(_dll, f'openmc_setting_get_{_suffix}')
    _getter.argtypes = [c_char_p, POINTER(_ctype)]
    _getter.restype = c_int
    _getter.errcheck = _error_handler

    _setter = getattr(_dll, f'openmc_setting_set_{_suffix}')
    _setter.argtypes = [c_char_p, _ctype]
    _setter.restype = c_int
    _setter.errcheck = _error_handler

    _SETTING_ACCESSORS[_ctype] = (_getter, _setter)

_dll.openmc_setting_get_string.argtypes = [c_char_p, POINTER(c_char_p)]
_dll.openmc_setting_get_string.restype = c_int
_dll.openmc_setting_get_string.errcheck = _error_handler


def _get_setting(ctype, name):
    value = ctype()
    getter, _ = _SETTING_ACCESSORS[ctype]
    getter(name.encode(), byref(value))
    return value.value


def _set_setting(ctype, name, value):
    _, setter = _SETTING_ACCESSORS[ctype]
    setter(name.encode(), value)


class _DLLFunctionProperty:
    """Data descriptor backed by type-specific C API setting functions."""

    def __init__(self, ctype, name):
        self.ctype = ctype
        self.name = name

    def __get__(self, instance, owner):
        if instance is None:
            return self
        return _get_setting(self.ctype, self.name)

    def __set__(self, instance, value):
        _set_setting(self.ctype, self.name, value)


class _Settings:
    # Attributes that are accessed through a descriptor
    cmfd_run = _DLLFunctionProperty(c_bool, 'cmfd_run')
    entropy_on = _DLLFunctionProperty(c_bool, 'entropy_on')
    event_based = _DLLFunctionProperty(c_bool, 'event_based')
    generations_per_batch = _DLLFunctionProperty(c_int32, 'gen_per_batch')
    inactive = _DLLFunctionProperty(c_int32, 'n_inactive')
    max_lost_particles = _DLLFunctionProperty(c_int32, 'max_lost_particles')
    max_write_lost_particles = _DLLFunctionProperty(
        c_int32, 'max_write_lost_particles'
    )
    need_depletion_rx = _DLLFunctionProperty(c_bool, 'need_depletion_rx')
    output_summary = _DLLFunctionProperty(c_bool, 'output_summary')
    particles = _DLLFunctionProperty(c_int64, 'n_particles')
    photon_transport = _DLLFunctionProperty(c_bool, 'photon_transport')
    rel_max_lost_particles = _DLLFunctionProperty(
        c_double, 'rel_max_lost_particles'
    )
    reduce_tallies = _DLLFunctionProperty(c_bool, 'reduce_tallies')
    restart_run = _DLLFunctionProperty(c_bool, 'restart_run')
    run_CE = _DLLFunctionProperty(c_bool, 'run_ce')
    trigger_on = _DLLFunctionProperty(c_bool, 'trigger_on')
    verbosity = _DLLFunctionProperty(c_int32, 'verbosity')
    weight_windows_on = _DLLFunctionProperty(c_bool, 'weight_windows_on')

    @property
    def run_mode(self):
        i = _get_setting(c_int32, 'run_mode')
        try:
            return _RUN_MODES[i]
        except KeyError:
            return None

    @run_mode.setter
    def run_mode(self, mode):
        for idx, mode_value in _RUN_MODES.items():
            if mode_value == mode:
                _set_setting(c_int32, 'run_mode', idx)
                break
        else:
            raise ValueError(f'Invalid run mode: {mode}')

    @property
    def path_statepoint(self):
        path = c_char_p()
        _dll.openmc_setting_get_string(b'path_statepoint', byref(path))
        return path.value.decode()

    @property
    def seed(self):
        return _dll.openmc_get_seed()

    @seed.setter
    def seed(self, seed):
        _dll.openmc_set_seed(seed)

    @property
    def stride(self):
        return _dll.openmc_get_stride()

    @stride.setter
    def stride(self, stride):
        _dll.openmc_set_stride(stride)

    def set_batches(self, n_batches, set_max_batches=True, add_sp_batch=True):
        """Set number of batches or maximum number of batches

        Parameters
        ----------
        n_batches : int
            Number of batches to simulate
        set_max_batches : bool
            Whether to set maximum number of batches. If True, the value of
            `n_max_batches` is overridden, otherwise the value of `n_batches`
            is overridden. Only has an effect when triggers are used
        add_sp_batch : bool
            Whether to add `n_batches` as a statepoint batch

        """
        _dll.openmc_set_n_batches(n_batches, set_max_batches, add_sp_batch)

    def get_batches(self, get_max_batches=True):
        """Get number of batches or maximum number of batches

        Parameters
        ----------
        get_max_batches : bool
            Return `n_max_batches` if true, else return `n_batches`. Difference
            arises only if triggers are used.

        Returns
        -------
        int
            Number of batches to simulate

        """
        n_batches = c_int()
        _dll.openmc_get_n_batches(n_batches, get_max_batches)

        return n_batches.value


settings = _Settings()
