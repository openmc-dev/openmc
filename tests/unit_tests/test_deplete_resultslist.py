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
    return openmc.deplete.ResultsList(filename)


def test_get_atoms(res):
    """Tests evaluating single nuclide concentration."""
    t, n = res.get_atoms("1", "Xe135")

    t_ref = [0.0, 1296000.0, 2592000.0, 3888000.0]
    n_ref = [6.674732823364922e+08, 4.022377919003244e+14,
             3.469325382312551e+14, 3.639463751297369e+14]

    np.testing.assert_allclose(t, t_ref)
    np.testing.assert_allclose(n, n_ref)


def test_get_reaction_rate(res):
    """Tests evaluating reaction rate."""
    t, r = res.get_reaction_rate("1", "Xe135", "(n,gamma)")

    t_ref = [0.0, 1296000.0, 2592000.0, 3888000.0]
    n_ref = [6.674732823364922e+08, 4.022377919003244e+14,
             3.469325382312551e+14, 3.639463751297369e+14]
    xs_ref = np.array([3.272630849911638e-05, 2.663333206501429e-05,
                       3.378999816491878e-05, 3.277013860196171e-05])

    np.testing.assert_allclose(t, t_ref)
    np.testing.assert_allclose(r, n_ref * xs_ref)


def test_get_eigenvalue(res):
    """Tests evaluating eigenvalue."""
    t, k = res.get_eigenvalue()

    t_ref = [0.0, 1296000.0, 2592000.0, 3888000.0]
    k_ref = [1.2162584685809341, 1.1021840248590875, 1.223040353281542,
             1.2360979537312164]

    np.testing.assert_allclose(t, t_ref)
    np.testing.assert_allclose(k, k_ref)
