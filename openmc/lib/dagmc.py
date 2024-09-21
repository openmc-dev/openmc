import sys
from ctypes import c_int, c_int32, POINTER, c_size_t

import numpy as np

from . import _dll
from .error import _error_handler


# DAGMC functions
_dll.openmc_get_dagmc_cell_ids.argtypes = [c_int32, POINTER(c_int32), POINTER(c_size_t)]
_dll.openmc_get_dagmc_cell_ids.restype = c_int
_dll.openmc_get_dagmc_cell_ids.errcheck = _error_handler
_dll.openmc_dagmc_universe_get_num_cells.argtypes = [c_int32, POINTER(c_size_t)]
_dll.openmc_dagmc_universe_get_num_cells.restype = c_int
_dll.openmc_dagmc_universe_get_num_cells.errcheck = _error_handler


def get_dagmc_cell_ids(dagmc_id):
    """Get the DAGMC cell IDs for a volume.

    Parameters
    ----------
    dagmc_id : int
        ID of the DAGMC Universe to get cell IDs from.
    n_cells : int
        Number of cells in the DAGMC Universe.

    Returns
    -------
    numpy.ndarray
        DAGMC cell IDs for the volume.

    """
    n = c_size_t()
    _dll.openmc_dagmc_universe_get_num_cells(dagmc_id, n)
    cell_ids = np.empty(n.value, dtype=np.int32)

    _dll.openmc_get_dagmc_cell_ids(
        dagmc_id,
        cell_ids.ctypes.data_as(POINTER(c_int32)),
        n
    )
    return cell_ids

def get_dagmc_universe_num_cells(dagmc_id):
    """Get the number of cells in a DAGMC universe.

    Parameters
    ----------
    dagmc_id : int
        ID of the DAGMC Universe to get the number of cell from.

    Returns
    -------
    int
        Number of cells in the DAGMC Universe.

    """
    n = c_size_t()
    _dll.openmc_dagmc_universe_get_num_cells(dagmc_id, n)
    return n.value
