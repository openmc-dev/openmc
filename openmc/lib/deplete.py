"""Ctypes bindings for C++ depletion solvers."""

from ctypes import c_int, c_double, POINTER

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
    c_int,          # solver_type
    _array_1d_dbl,  # result
]


def cram_solve(matrix, n0, dt, order=48, solver_type=0):
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
    solver_type : int
        Solver type: 0 for general LU (only type currently supported).

    Returns
    -------
    numpy.ndarray
        Final atom numbers.

    """
    from scipy.sparse import csc_array
    if order not in (16, 48):
        raise ValueError(f"CRAM order must be 16 or 48, got {order}")
    if solver_type not in (0,):
        raise ValueError(
            f"Only solver_type 0 (general LU) is supported, got {solver_type}")

    matrix = csc_array(matrix)
    n = matrix.shape[0]
    indptr = np.asarray(matrix.indptr, dtype=np.int32)
    indices = np.asarray(matrix.indices, dtype=np.int32)
    data = np.asarray(matrix.data, dtype=np.float64)
    n0 = np.asarray(n0, dtype=np.float64)
    result = np.empty(n, dtype=np.float64)

    _dll.openmc_cram_solve(
        n, indptr, indices, data, n0, dt, order, solver_type, result)

    return result


# --- CRAM batch solve ---

_dll.openmc_cram_solve_batch.restype = c_int
_dll.openmc_cram_solve_batch.errcheck = _error_handler
_dll.openmc_cram_solve_batch.argtypes = [
    c_int,          # n_materials
    _array_1d_int,  # dims
    _array_1d_int,  # all_indptr
    _array_1d_int,  # all_indices
    _array_1d_dbl,  # all_data
    _array_1d_int,  # nnz_per_mat
    _array_1d_dbl,  # all_n0
    c_double,       # dt
    c_int,          # order
    c_int,          # solver_type
    _array_1d_dbl,  # all_results
]


def cram_solve_batch(matrices, n0_list, dt, order=48, solver_type=0):
    """Solve multiple Bateman systems in parallel using C++ CRAM with OpenMP.

    Each system is defined by a sparse transmutation matrix and an initial
    atom number vector. All systems are solved over the same time step.
    The matrices may have different dimensions and sparsity patterns.

    Parameters
    ----------
    matrices : list of scipy.sparse.csc_array
        Sparse transmutation matrices, one per material.
    n0_list : list of numpy.ndarray
        Initial atom number vectors, one per material. ``n0_list[m]`` must
        have length equal to ``matrices[m].shape[0]``.
    dt : float
        Time step in seconds (shared by all materials).
    order : int
        CRAM approximation order (16 or 48).
    solver_type : int
        Solver type: 0 for general LU (only type currently supported).

    Returns
    -------
    list of numpy.ndarray
        Final atom numbers for each material.

    """
    if order not in (16, 48):
        raise ValueError(f"CRAM order must be 16 or 48, got {order}")
    if solver_type not in (0,):
        raise ValueError(
            f"Only solver_type 0 (general LU) is supported, got {solver_type}")

    n_materials = len(matrices)
    dims = np.array([m.shape[0] for m in matrices], dtype=np.int32)
    nnz_per_mat = np.array(
        [m.indptr[-1] for m in matrices], dtype=np.int32)

    all_indptr = np.concatenate(
        [np.asarray(m.indptr, dtype=np.int32) for m in matrices])
    all_indices = np.concatenate(
        [np.asarray(m.indices, dtype=np.int32) for m in matrices])
    all_data = np.concatenate(
        [np.asarray(m.data, dtype=np.float64) for m in matrices])
    all_n0 = np.concatenate(
        [np.asarray(v, dtype=np.float64) for v in n0_list])
    all_results = np.empty(dims.sum(), dtype=np.float64)

    _dll.openmc_cram_solve_batch(
        n_materials, dims, all_indptr, all_indices, all_data,
        nnz_per_mat, all_n0, dt, order, solver_type, all_results)

    # Split concatenated results back into per-material arrays
    offsets = np.cumsum(dims)
    return [all_results[s:e] for s, e in zip(
        np.concatenate([[0], offsets[:-1]]), offsets)]
