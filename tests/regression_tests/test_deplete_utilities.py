""" Tests the utilities classes.

This also tests the results read/write code.
"""

from pathlib import Path

import numpy as np
import pytest
from openmc.deplete import results
from openmc.deplete import utilities


@pytest.fixture
def res():
    """Load the reference results"""
    filename = str(Path(__file__).with_name('test_reference.h5'))
    return results.read_results(filename)


def test_evaluate_single_nuclide(res):
    """Tests evaluating single nuclide utility code."""
    x, y = utilities.evaluate_single_nuclide(res, "1", "Xe135")

    x_ref = [0.0, 1296000.0, 2592000.0, 3888000.0]
    y_ref = [6.6747328233649218e+08, 3.5519299354458244e+14,
             3.4599104054580338e+14, 3.3821165110278112e+14]

    np.testing.assert_array_equal(x, x_ref)
    np.testing.assert_array_equal(y, y_ref)

def test_evaluate_reaction_rate(res):
    """Tests evaluating reaction rate utility code."""
    x, y = utilities.evaluate_reaction_rate(res, "1", "Xe135", "(n,gamma)")

    x_ref = [0.0, 1296000.0, 2592000.0, 3888000.0]
    xe_ref = np.array([6.6747328233649218e+08, 3.5519299354458244e+14,
                       3.4599104054580338e+14, 3.3821165110278112e+14])
    r_ref = np.array([4.0643598574337784e-05, 4.1457730544386974e-05,
                      3.4121248544056681e-05, 3.9204686657643301e-05])

    np.testing.assert_array_equal(x, x_ref)
    np.testing.assert_array_equal(y, xe_ref * r_ref)


def test_evaluate_eigenvalue(res):
    """Tests evaluating eigenvalue."""
    x, y = utilities.evaluate_eigenvalue(res)

    x_ref = [0.0, 1296000.0, 2592000.0, 3888000.0]
    y_ref = [1.1921986054449838, 1.1712785643938586, 1.1927099024502694, 1.2269183590698847]

    np.testing.assert_array_equal(x, x_ref)
    np.testing.assert_array_equal(y, y_ref)
