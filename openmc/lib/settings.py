from ctypes import (c_int, c_int32, c_int64, c_double, c_char_p, c_bool,
                    POINTER)

from . import _dll
from .core import _DLLGlobal
from .error import _error_handler

_RUN_MODES = {0: 'fixed source',
              1: 'eigenvalue',
              2: 'plot',
              3: 'particle restart',
              4: 'volume'}

_dll.openmc_set_seed.argtypes = [c_int64]
_dll.openmc_get_seed.restype = c_int64


class _Settings(object):
    # Attributes that are accessed through a descriptor
    batches = _DLLGlobal(c_int32, 'n_batches')
    cmfd_run = _DLLGlobal(c_bool, 'cmfd_run')
    entropy_on = _DLLGlobal(c_bool, 'entropy_on')
    generations_per_batch = _DLLGlobal(c_int32, 'gen_per_batch')
    inactive = _DLLGlobal(c_int32, 'n_inactive')
    particles = _DLLGlobal(c_int64, 'n_particles')
    restart_run = _DLLGlobal(c_bool, 'restart_run')
    run_CE = _DLLGlobal(c_bool, 'run_CE')
    verbosity = _DLLGlobal(c_int, 'verbosity')
    output_summary = _DLLGlobal(c_bool, 'output_summary')

    @property
    def run_mode(self):
        raise Exception("Setting run_mode through python API is no longer supported.")

    @run_mode.setter
    def run_mode(self, mode):
        raise Exception("Setting run_mode through python API is no longer supported.")

    @property
    def path_statepoint(self):
        path = c_char_p.in_dll(_dll, 'path_statepoint').value
        return path.decode()

    @property
    def seed(self):
        return _dll.openmc_get_seed()

    @seed.setter
    def seed(self, seed):
        _dll.openmc_set_seed(seed)


settings = _Settings()
