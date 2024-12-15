from ctypes import c_int, c_int32, c_int64, c_double, c_char_p, c_bool, POINTER

from . import _dll
from .core import _DLLGlobal
from .error import _error_handler

_RUN_MODES = {
    1: "fixed source",
    2: "eigenvalue",
    3: "plot",
    4: "particle restart",
    5: "volume",
}

_dll.openmc_set_seed.argtypes = [c_int64]
_dll.openmc_get_seed.restype = c_int64
_dll.openmc_get_n_batches.argtypes = [POINTER(c_int), c_bool]
_dll.openmc_get_n_batches.restype = c_int
_dll.openmc_get_n_batches.errcheck = _error_handler
_dll.openmc_set_n_batches.argtypes = [c_int32, c_bool, c_bool]
_dll.openmc_set_n_batches.restype = c_int
_dll.openmc_set_n_batches.errcheck = _error_handler


class _Settings:
    # Attributes that are accessed through a descriptor
    cmfd_run = _DLLGlobal(c_bool, "cmfd_run")
    entropy_on = _DLLGlobal(c_bool, "entropy_on")
    generations_per_batch = _DLLGlobal(c_int32, "gen_per_batch")
    inactive = _DLLGlobal(c_int32, "n_inactive")
    max_lost_particles = _DLLGlobal(c_int32, "max_lost_particles")
    need_depletion_rx = _DLLGlobal(c_bool, "need_depletion_rx")
    output_summary = _DLLGlobal(c_bool, "output_summary")
    particles = _DLLGlobal(c_int64, "n_particles")
    rel_max_lost_particles = _DLLGlobal(c_double, "rel_max_lost_particles")
    restart_run = _DLLGlobal(c_bool, "restart_run")
    run_CE = _DLLGlobal(c_bool, "run_CE")
    verbosity = _DLLGlobal(c_int, "verbosity")
    event_based = _DLLGlobal(c_bool, "event_based")
    weight_windows_on = _DLLGlobal(c_bool, "weight_windows_on")

    @property
    def run_mode(self):
        i = c_int.in_dll(_dll, "run_mode").value
        try:
            return _RUN_MODES[i]
        except KeyError:
            return None

    @run_mode.setter
    def run_mode(self, mode):
        current_idx = c_int.in_dll(_dll, "run_mode")
        for idx, mode_value in _RUN_MODES.items():
            if mode_value == mode:
                current_idx.value = idx
                break
        else:
            raise ValueError(f"Invalid run mode: {mode}")

    @property
    def path_statepoint(self):
        path = c_char_p.in_dll(_dll, "path_statepoint_c").value
        return path.decode()

    @property
    def seed(self):
        return _dll.openmc_get_seed()

    @seed.setter
    def seed(self, seed):
        _dll.openmc_set_seed(seed)

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
