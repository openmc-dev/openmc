"""Tests the ResultsList class"""

from pathlib import Path

import numpy as np
import pytest
import openmc.deplete


@pytest.fixture
def res():
    """Load the reference results"""
    filename = Path(__file__).parents[1] / 'regression_tests' / 'test_reference.h5'
    return openmc.deplete.ResultsList(filename)


def test_get_atoms(res):
    """Tests evaluating single nuclide concentration."""
    t, n = res.get_atoms("1", "Xe135")

    t_ref = [0.0, 1296000.0, 2592000.0, 3888000.0]
    n_ref = [6.6747328233649218e+08, 3.4589992012016338e+14,
             3.5635060369969225e+14, 3.5195113630100188e+14]

    np.testing.assert_allclose(t, t_ref)
    np.testing.assert_allclose(n, n_ref)


def test_get_reaction_rate(res):
    """Tests evaluating reaction rate."""
    t, r = res.get_reaction_rate("1", "Xe135", "(n,gamma)")

    t_ref = [0.0, 1296000.0, 2592000.0, 3888000.0]
    n_ref = np.array([6.6747328233649218e+08, 3.4589992012016338e+14,
                      3.5635060369969225e+14, 3.5195113630100188e+14])
    xs_ref = np.array([4.1340608491478010e-05, 4.1120938620476115e-05,
                       4.3341529708654921e-05, 3.8716623651147821e-05])

    np.testing.assert_allclose(t, t_ref)
    np.testing.assert_allclose(r, n_ref * xs_ref)


def test_get_eigenvalue(res):
    """Tests evaluating eigenvalue."""
    t, k = res.get_eigenvalue()

    t_ref = [0.0, 1296000.0, 2592000.0, 3888000.0]
    k_ref = [1.1798617938070866, 1.1745713141097096, 1.1732427763487678,
             1.213699703239334]

    np.testing.assert_allclose(t, t_ref)
    np.testing.assert_allclose(k, k_ref)
