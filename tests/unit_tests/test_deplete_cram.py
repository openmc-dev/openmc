"""Tests for cram.py.

Compares a few Mathematica matrix exponentials to CRAM16/CRAM48.
Tests substep accuracy against self-converged reference solutions.
"""

import numpy as np
import pytest
import scipy.sparse as sp
from pytest import approx

from openmc.deplete.cram import (CRAM16, CRAM48, Cram16Solver, Cram48Solver,
                                  IPFCramSolver)
from openmc.lib.deplete import cram_solve


def test_CRAM16():
    """Test 16th order CRAM against Mathematica reference."""
    x = np.array([1.0, 1.0])
    mat = sp.csc_array([[-1.0, 0.0], [-2.0, -3.0]])
    dt = 0.1

    z = CRAM16(mat, x, dt)

    z0 = np.array((0.904837418035960, 0.576799023327476))

    assert z == approx(z0)


def test_CRAM48():
    """Test 48th order CRAM against Mathematica reference."""
    x = np.array([1.0, 1.0])
    mat = sp.csc_array([[-1.0, 0.0], [-2.0, -3.0]])
    dt = 0.1

    z = CRAM48(mat, x, dt)

    z0 = np.array((0.904837418035960, 0.576799023327476))

    assert z == approx(z0)


def test_cram_solver_objects():
    """Develop API singleton solver objects remain available and callable."""
    assert isinstance(Cram16Solver, IPFCramSolver)
    assert isinstance(Cram48Solver, IPFCramSolver)

    x = np.array([1.0, 1.0])
    mat = sp.csc_array([[-1.0, 0.0], [-2.0, -3.0]])
    dt = 0.1

    assert Cram16Solver(mat, x, dt) == approx(CRAM16(mat, x, dt))
    assert Cram48Solver(mat, x, dt) == approx(CRAM48(mat, x, dt))


def test_cram_solve_cpp():
    """Test C++ cram_solve binding directly against Mathematica reference."""
    A = sp.csc_array([[-1.0, 0.0], [-2.0, -3.0]])
    n0 = np.array([1.0, 1.0])
    dt = 0.1

    z48 = cram_solve(A, n0, dt, order=48)
    z16 = cram_solve(A, n0, dt, order=16)

    z_ref = np.array((0.904837418035960, 0.576799023327476))

    assert z48 == approx(z_ref)
    assert z16 == approx(z_ref)


def test_substeps1_matches_original():
    """substeps=1 must be bitwise identical to original result."""
    x = np.array([1.0, 1.0])
    mat = sp.csc_array([[-1.0, 0.0], [-2.0, -3.0]])
    dt = 0.1

    z_orig = CRAM48(mat, x, dt)
    z_sub1 = CRAM48(mat, x, dt, substeps=1)
    z_cpp = cram_solve(mat, x, dt, order=48, substeps=1)

    np.testing.assert_array_equal(z_sub1, z_orig)
    np.testing.assert_array_equal(z_cpp, z_orig)


def test_substeps2_matches_two_half_steps():
    """substeps=2 must match two independent CRAM calls with dt/2."""
    x = np.array([1.0, 1.0])
    mat = sp.csc_array([[-1.0, 0.0], [-2.0, -3.0]])
    dt = 1.0

    z_half = CRAM48(mat, x, dt / 2)
    z_two = CRAM48(mat, z_half, dt / 2)
    z_sub2 = CRAM48(mat, x, dt, substeps=2)
    z_cpp = cram_solve(mat, x, dt, order=48, substeps=2)

    assert z_sub2 == approx(z_two, rel=1e-12)
    assert z_cpp == approx(z_two, rel=1e-12)


@pytest.mark.parametrize("substeps", [0, -1])
def test_invalid_substeps(substeps):
    """substeps must be a positive integer at call time."""
    x = np.array([1.0, 1.0])
    mat = sp.csc_array([[-1.0, 0.0], [-2.0, -3.0]])

    with pytest.raises(ValueError, match="substeps"):
        CRAM48(mat, x, 0.1, substeps=substeps)

    with pytest.raises(ValueError, match="substeps"):
        cram_solve(mat, x, 0.1, order=48, substeps=substeps)


def test_substeps_self_convergence():
    """Increasing substeps converges toward reference solution.

    Uses CRAM16 where substep convergence is visible.
    """
    mat = sp.csc_array([[-1.0, 0.0], [-2.0, -3.0]])
    x = np.array([1.0, 1.0])
    dt = 50.0

    n_ref = CRAM16(mat, x, dt, substeps=128)

    prev_err = np.inf
    for substeps in [1, 2, 4, 8, 16]:
        n_sub = CRAM16(mat, x, dt, substeps=substeps)
        err = np.linalg.norm(n_sub - n_ref) / np.linalg.norm(n_ref)
        assert err < prev_err
        prev_err = err
