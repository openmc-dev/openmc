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


def test_CRAM16():
    """Test 16-term CRAM."""
    x = np.array([1.0, 1.0])
    mat = sp.csr_matrix([[-1.0, 0.0], [-2.0, -3.0]])
    dt = 0.1

    z = CRAM16(mat, x, dt)

    # Solution from mathematica
    z0 = np.array((0.904837418035960, 0.576799023327476))

    assert z == approx(z0)


def test_CRAM48():
    """Test 48-term CRAM."""
    x = np.array([1.0, 1.0])
    mat = sp.csr_matrix([[-1.0, 0.0], [-2.0, -3.0]])
    dt = 0.1

    z = CRAM48(mat, x, dt)

    # Solution from mathematica
    z0 = np.array((0.904837418035960, 0.576799023327476))

    assert z == approx(z0)


def test_substeps1_matches_original():
    """substeps=1 must be bitwise identical to original spsolve path."""
    x = np.array([1.0, 1.0])
    mat = sp.csr_matrix([[-1.0, 0.0], [-2.0, -3.0]])
    dt = 0.1

    z_orig = CRAM48(mat, x, dt)
    z_sub1 = CRAM48(mat, x, dt, substeps=1)

    np.testing.assert_array_equal(z_sub1, z_orig)


def test_substeps2_matches_two_half_steps():
    """substeps=2 must match two independent CRAM calls with dt/2."""
    x = np.array([1.0, 1.0])
    mat = sp.csr_matrix([[-1.0, 0.0], [-2.0, -3.0]])
    dt = 1.0

    # Two manual half-steps using original spsolve path
    z_half = CRAM48(mat, x, dt / 2)
    z_two = CRAM48(mat, z_half, dt / 2)

    # Single call with substeps=2
    z_sub2 = CRAM48(mat, x, dt, substeps=2)

    assert z_sub2 == approx(z_two, rel=1e-12)


@pytest.mark.parametrize("substeps", [0, -1])
def test_invalid_substeps(substeps):
    """substeps must be a positive integer at call time."""
    x = np.array([1.0, 1.0])
    mat = sp.csr_matrix([[-1.0, 0.0], [-2.0, -3.0]])

    with pytest.raises(ValueError, match="substeps"):
        CRAM48(mat, x, 0.1, substeps=substeps)


def test_substeps_self_convergence():
    """Increasing substeps converges toward reference solution.

    Uses CRAM16 (alpha0 ~ 2e-16) where substep convergence is visible.
    CRAM48 (alpha0 ~ 2e-47) is already near machine precision for small
    systems; its correctness is verified by the other substep tests.
    """
    mat = sp.csr_matrix([[-1.0, 0.0], [-2.0, -3.0]])
    x = np.array([1.0, 1.0])
    dt = 50  # lambda*dt = 50 and 150, stresses CRAM16

    n_ref = CRAM16(mat, x, dt, substeps=128)

    prev_err = np.inf
    for s in [1, 2, 4, 8, 16]:
        n_s = CRAM16(mat, x, dt, substeps=s)
        err = np.linalg.norm(n_s - n_ref) / np.linalg.norm(n_ref)
        assert err < prev_err, \
            f"substeps={s} error {err:.2e} not less than previous {prev_err:.2e}"
        prev_err = err
