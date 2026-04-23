"""Chebyshev Rational Approximation Method module

Provides CRAM16 and CRAM48 solvers for use in openmc.deplete, backed by the
C++ IPFCramSolver implementation.
"""

import numpy as np

from .._sparse_compat import csc_array
from openmc.lib.deplete import cram_solve

__all__ = ["CRAM16", "CRAM48"]


def _cram_solve(A, n0, dt, order, substeps=1):
    """Single-material CRAM solve via C++ IPFCramSolver."""
    return cram_solve(
        csc_array(A, dtype=np.float64), np.asarray(n0, dtype=np.float64),
        float(dt), order=order, substeps=substeps)


def CRAM48(A, n0, dt, substeps=1):
    """Solve depletion equations using 48th order IPF CRAM

    Parameters
    ----------
    A : scipy.sparse.csc_array
        Sparse transmutation matrix ``A[j, i]`` describing rates at
        which isotope ``i`` transmutes to isotope ``j``
    n0 : numpy.ndarray
        Initial compositions, typically given in number of atoms in some
        material or an atom density
    dt : float
        Time [s] of the specific interval to be solved
    substeps : int, optional
        Number of equal substeps to use within ``dt``

    Returns
    -------
    numpy.ndarray
        Final compositions after ``dt``

    """
    return _cram_solve(A, n0, dt, order=48, substeps=substeps)


def CRAM16(A, n0, dt, substeps=1):
    """Solve depletion equations using 16th order IPF CRAM

    Parameters
    ----------
    A : scipy.sparse.csc_array
        Sparse transmutation matrix ``A[j, i]`` describing rates at
        which isotope ``i`` transmutes to isotope ``j``
    n0 : numpy.ndarray
        Initial compositions, typically given in number of atoms in some
        material or an atom density
    dt : float
        Time [s] of the specific interval to be solved
    substeps : int, optional
        Number of equal substeps to use within ``dt``

    Returns
    -------
    numpy.ndarray
        Final compositions after ``dt``

    """
    return _cram_solve(A, n0, dt, order=16, substeps=substeps)
