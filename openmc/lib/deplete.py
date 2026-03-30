"""Ctypes bindings for C++ depletion solvers."""

from ctypes import c_int, c_double

import numpy as np
from numpy.ctypeslib import ndpointer

from .error import _error_handler
from . import _dll

_array_1d_int = ndpointer(dtype=np.int32, ndim=1, flags='CONTIGUOUS')
_array_1d_dbl = ndpointer(dtype=np.float64, ndim=1, flags='CONTIGUOUS')

# --- CRAM single-material solve ---

_dll.openmc_cram_solve.restype = c_int
_dll.openmc_cram_solve.errcheck = _error_handler
_dll.openmc_cram_solve.argtypes = [
    c_int,          # n
    _array_1d_int,  # indptr
    _array_1d_int,  # indices
    _array_1d_dbl,  # data
    _array_1d_dbl,  # n0
    c_double,       # dt
    c_int,          # order
    _array_1d_dbl,  # result
]


def cram_solve(matrix, n0, dt, order=48):
    """Solve a single Bateman system using C++ CRAM.

    Parameters
    ----------
    matrix : scipy.sparse.csc_array
        Sparse transmutation matrix.
    n0 : numpy.ndarray
        Initial atom number vector.
    dt : float
        Time step in seconds.
    order : int
        CRAM approximation order (16 or 48).

    Returns
    -------
    numpy.ndarray
        Final atom numbers.

    """
    from scipy.sparse import csc_array
    if order not in (16, 48):
        raise ValueError(f"CRAM order must be 16 or 48, got {order}")

    matrix = csc_array(matrix)
    n = matrix.shape[0]
    indptr = np.asarray(matrix.indptr, dtype=np.int32)
    indices = np.asarray(matrix.indices, dtype=np.int32)
    data = np.asarray(matrix.data, dtype=np.float64)
    n0 = np.asarray(n0, dtype=np.float64)
    result = np.empty(n, dtype=np.float64)

    _dll.openmc_cram_solve(
        n, indptr, indices, data, n0, dt, order, result)

    return result
