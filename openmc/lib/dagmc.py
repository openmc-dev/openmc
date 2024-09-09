import sys
from ctypes import c_int, c_int32, POINTER, c_size_t

import numpy as np

from . import _dll
from .error import _error_handler


# DAGMC functions
_dll.openmc_get_dagmc_cell_ids.argtypes = [c_int32, POINTER(c_int32), POINTER(c_size_t)]
_dll.openmc_get_dagmc_cell_ids.restype = c_int
_dll.openmc_get_dagmc_cell_ids.errcheck = _error_handler


def get_dagmc_cell_ids(volume_id, n_cells):
    """get the dagmc cell ids for a volume.

    parameters
    ----------
    volume_id : int
        id of the volume to get dagmc cell ids for.
    n_cells : int
        number of cells in the volume.

    returns
    -------
    numpy.ndarray
        dagmc cell ids for the volume.

    """
    cell_ids = np.empty(n_cells, dtype=np.int32)
    n = c_size_t()
    _dll.openmc_get_dagmc_cell_ids(
        volume_id,
        cell_ids.ctypes.data_as(POINTER(c_int32)),
        n
    )
    if n.value != n_cells:
        raise ValueError(f"Number of cells obtained {n.value} from dagmc does "
                         f"not match the expected number of cells {n_cells}.")
    return cell_ids