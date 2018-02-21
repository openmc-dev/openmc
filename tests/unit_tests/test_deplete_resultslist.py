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
    n_ref = [6.6747328233649218e+08, 3.5421791038348462e+14,
             3.6208592242443462e+14, 3.3799758969347038e+14]

    np.testing.assert_array_equal(t, t_ref)
    np.testing.assert_array_equal(n, n_ref)

def test_get_reaction_rate(res):
    """Tests evaluating reaction rate."""
    t, r = res.get_reaction_rate("1", "Xe135", "(n,gamma)")

    t_ref = [0.0, 1296000.0, 2592000.0, 3888000.0]
    n_ref = np.array([6.6747328233649218e+08, 3.5421791038348462e+14,
                       3.6208592242443462e+14, 3.3799758969347038e+14])
    xs_ref = np.array([4.0594392323131994e-05, 3.9249546927524987e-05,
                       3.8394587728581798e-05, 4.1521845978371697e-05])

    np.testing.assert_array_equal(t, t_ref)
    np.testing.assert_array_equal(r, n_ref * xs_ref)


def test_get_eigenvalue(res):
    """Tests evaluating eigenvalue."""
    t, k = res.get_eigenvalue()

    t_ref = [0.0, 1296000.0, 2592000.0, 3888000.0]
    k_ref = [1.181281798790367, 1.1798750921988739, 1.1965943696058159,
             1.2207119847790813]

    np.testing.assert_array_equal(t, t_ref)
    np.testing.assert_array_equal(k, k_ref)
