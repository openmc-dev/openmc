"""Chebyshev Rational Approximation Method module.

Provides :func:`CRAM16` and :func:`CRAM48` solvers for use in
:mod:`openmc.deplete`, backed by the C++ ``IPFCramSolver`` implementation in
:mod:`openmc.lib.deplete`.
"""

import numpy as np

from openmc.lib.deplete import cram_solve
from .._sparse_compat import csc_array

__all__ = ["CRAM16", "CRAM48"]


def _cram_solve(A, n0, dt, order, substeps=1):
    """Single-material CRAM solve via the C++ ``IPFCramSolver`` backend."""
    return cram_solve(
        csc_array(A, dtype=np.float64), np.asarray(n0, dtype=np.float64),
        float(dt), order=order, substeps=substeps)


def CRAM48(A, n0, dt, substeps=1):
    r"""Solve depletion equations using 48th order IPF CRAM.

    Implements the Incomplete Partial Fraction form of the Chebyshev
    Rational Approximation Method (CRAM) as described in M. Pusa,
    "`Higher-Order Chebyshev Rational Approximation Method and Application to
    Burnup Equations <https://doi.org/10.13182/NSE15-26>`_," Nucl. Sci. Eng.,
    182:3, 297-318. When ``substeps > 1`` the time interval is split into
    equal sub-intervals and LU factorizations are reused across them as
    described in A. Isotalo and M. Pusa, "`Improving the Accuracy of the
    Chebyshev Rational Approximation Method Using Substeps
    <https://doi.org/10.13182/NSE15-67>`_," Nucl. Sci. Eng., 183:1, 65-77.

    Parameters
    ----------
    A : scipy.sparse.csc_array
        Sparse transmutation matrix ``A[j, i]`` describing rates at
        which isotope ``i`` transmutes to isotope ``j``.
    n0 : numpy.ndarray
        Initial compositions, typically given in number of atoms in some
        material or an atom density.
    dt : float
        Time [s] of the interval to be solved.
    substeps : int, optional
        Number of equal substeps to take within ``dt``. Defaults to 1.

    Returns
    -------
    numpy.ndarray
        Final compositions after ``dt``.

    """
    return _cram_solve(A, n0, dt, order=48, substeps=substeps)


def CRAM16(A, n0, dt, substeps=1):
    r"""Solve depletion equations using 16th order IPF CRAM.

    See :func:`CRAM48` for a description of the method and substep behavior.

    Parameters
    ----------
    A : scipy.sparse.csc_array
        Sparse transmutation matrix ``A[j, i]`` describing rates at
        which isotope ``i`` transmutes to isotope ``j``.
    n0 : numpy.ndarray
        Initial compositions, typically given in number of atoms in some
        material or an atom density.
    dt : float
        Time [s] of the interval to be solved.
    substeps : int, optional
        Number of equal substeps to take within ``dt``. Defaults to 1.

    Returns
    -------
    numpy.ndarray
        Final compositions after ``dt``.

    """
    return _cram_solve(A, n0, dt, order=16, substeps=substeps)
