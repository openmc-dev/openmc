"""Tests the ResultsList class"""

from pathlib import Path

import numpy as np
import pytest
import openmc.deplete


@pytest.fixture
def res():
    """Load the reference results"""
    filename = (Path(__file__).parents[1] / 'regression_tests' / 'deplete'
                / 'test_reference.h5')
    return openmc.deplete.ResultsList.from_hdf5(filename)


def test_get_atoms(res):
    """Tests evaluating single nuclide concentration."""
    t, n = res.get_atoms("1", "Xe135")

    t_ref = [0.0, 1296000.0, 2592000.0, 3888000.0]
    n_ref = [6.67473282e+08, 3.76986925e+14, 3.68587383e+14, 3.91338675e+14]

    np.testing.assert_allclose(t, t_ref)
    np.testing.assert_allclose(n, n_ref)


def test_get_reaction_rate(res):
    """Tests evaluating reaction rate."""
    t, r = res.get_reaction_rate("1", "Xe135", "(n,gamma)")

    t_ref = [0.0, 1296000.0, 2592000.0, 3888000.0]
    n_ref = [6.67473282e+08, 3.76986925e+14, 3.68587383e+14, 3.91338675e+14]
    xs_ref = [3.32282266e-05, 2.76207120e-05, 4.10986677e-05, 3.72453665e-05]

    np.testing.assert_allclose(t, t_ref)
    np.testing.assert_allclose(r, np.array(n_ref) * xs_ref)


def test_get_eigenvalue(res):
    """Tests evaluating eigenvalue."""
    t, k = res.get_eigenvalue()

    t_ref = [0.0, 1296000.0, 2592000.0, 3888000.0]
    k_ref = [1.16984322, 1.19097427, 1.03012572, 1.20045627]
    u_ref = [0.0375587, 0.0347639, 0.07216021, 0.02839642]

    np.testing.assert_allclose(t, t_ref)
    np.testing.assert_allclose(k[:, 0], k_ref)
    np.testing.assert_allclose(k[:, 1], u_ref)
