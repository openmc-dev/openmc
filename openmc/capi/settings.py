from ctypes import c_int, c_int32, c_int64, c_double, c_char_p, POINTER

from . import _dll
from .error import _error_handler

_RUN_MODES = {1: 'fixed source',
              2: 'eigenvalue',
              3: 'plot',
              4: 'particle restart',
              5: 'volume'}


class _Settings(object):
    @property
    def batches(self):
        return c_int32.in_dll(_dll, 'n_batches').value

    @batches.setter
    def batches(self, n):
        _dll.openmc_set_batches(n)

    @property
    def generations_per_batch(self):
        return c_int32.in_dll(_dll, 'gen_per_batch').value

    @generations_per_batch.setter
    def generations_per_batch(self, n):
        _dll.openmc_set_generations_per_batch(n)

    @property
    def inactive(self):
        return c_int32.in_dll(_dll, 'n_inactive').value

    @inactive.setter
    def inactive(self, n):
        _dll.openmc_set_inactive_batches(n)

    @property
    def particles(self):
        return c_int64.in_dll(_dll, 'n_particles').value

    @particles.setter
    def particles(self, n):
        _dll.openmc_set_particles(n)

    @property
    def run_mode(self):
        i = c_int.in_dll(_dll, 'run_mode').value
        try:
            return _RUN_MODES[i]
        except KeyError:
            return None

    @run_mode.setter
    def run_mode(self, mode):
        current_idx = c_int.in_dll(_dll, 'run_mode')
        for idx, mode_value in _RUN_MODES.items():
            if mode_value == mode:
                current_idx.value = idx
                break
        else:
            raise ValueError('Invalid run mode: {}'.format(mode))


settings = _Settings()
