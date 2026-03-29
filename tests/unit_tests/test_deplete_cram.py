""" Tests for cram.py

Compares a few Mathematica matrix exponentials to CRAM16/CRAM48.
Tests substep accuracy against self-converged reference solutions.
"""

from pytest import approx
import numpy as np
import scipy.sparse as sp
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


# --- Substep tests ---

def test_substeps_default():
    """Default substeps=1 on module-level solvers."""
    assert Cram48Solver.substeps == 1
    assert Cram16Solver.substeps == 1


def test_substeps1_matches_original():
    """substeps=1 must be bitwise identical to original spsolve path."""
    x = np.array([1.0, 1.0])
    mat = sp.csr_matrix([[-1.0, 0.0], [-2.0, -3.0]])
    dt = 0.1

    z_orig = CRAM48(mat, x, dt)
    solver_1 = IPFCramSolver(Cram48Solver.alpha, Cram48Solver.theta,
                              Cram48Solver.alpha0, substeps=1)
    z_sub1 = solver_1(mat, x, dt)

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
    solver_2 = IPFCramSolver(Cram48Solver.alpha, Cram48Solver.theta,
                              Cram48Solver.alpha0, substeps=2)
    z_sub2 = solver_2(mat, x, dt)

    assert z_sub2 == approx(z_two, rel=1e-12)


def test_substeps_self_convergence():
    """Increasing substeps converges toward reference solution.

    Uses CRAM16 (alpha0 ~ 2e-16) where substep convergence is visible.
    CRAM48 (alpha0 ~ 2e-47) is already near machine precision for small
    systems; its correctness is verified by the other substep tests.
    """
    mat = sp.csr_matrix([[-1.0, 0.0], [-2.0, -3.0]])
    x = np.array([1.0, 1.0])
    dt = 50  # lambda*dt = 50 and 150, stresses CRAM16

    ref_solver = IPFCramSolver(Cram16Solver.alpha, Cram16Solver.theta,
                                Cram16Solver.alpha0, substeps=128)
    n_ref = ref_solver(mat, x, dt)

    prev_err = np.inf
    for s in [1, 2, 4, 8, 16]:
        solver = IPFCramSolver(Cram16Solver.alpha, Cram16Solver.theta,
                                Cram16Solver.alpha0, substeps=s)
        n_s = solver(mat, x, dt)
        err = np.linalg.norm(n_s - n_ref) / np.linalg.norm(n_ref)
        assert err < prev_err, \
            f"substeps={s} error {err:.2e} not less than previous {prev_err:.2e}"
        prev_err = err
