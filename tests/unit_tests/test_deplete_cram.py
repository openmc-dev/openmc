"""Tests for cram.py

Compares Mathematica matrix exponentials to CRAM16/CRAM48 (C++ backend)
and tests the C++ cram_solve binding directly.
"""

from pytest import approx
import numpy as np
import scipy.sparse as sp
from openmc.deplete.cram import CRAM16, CRAM48


def test_CRAM16():
    """Test 16th order CRAM against Mathematica reference."""
    x = np.array([1.0, 1.0])
    mat = sp.csc_array([[-1.0, 0.0], [-2.0, -3.0]])
    dt = 0.1

    z = CRAM16(mat, x, dt)

    # Solution from Mathematica
    z0 = np.array((0.904837418035960, 0.576799023327476))

    assert z == approx(z0)


def test_CRAM48():
    """Test 48th order CRAM against Mathematica reference."""
    x = np.array([1.0, 1.0])
    mat = sp.csc_array([[-1.0, 0.0], [-2.0, -3.0]])
    dt = 0.1

    z = CRAM48(mat, x, dt)

    # Solution from Mathematica
    z0 = np.array((0.904837418035960, 0.576799023327476))

    assert z == approx(z0)


def test_cram_solve_cpp():
    """Test C++ cram_solve binding directly against Mathematica reference."""
    from openmc.lib.deplete import cram_solve

    A = sp.csc_array([[-1.0, 0.0], [-2.0, -3.0]])
    n0 = np.array([1.0, 1.0])
    dt = 0.1

    z48 = cram_solve(A, n0, dt, order=48)
    z16 = cram_solve(A, n0, dt, order=16)

    # Reference from Mathematica
    z_ref = np.array((0.904837418035960, 0.576799023327476))

    assert z48 == approx(z_ref)
    assert z16 == approx(z_ref)
