""" Tests for cram.py

Compares a few Mathematica matrix exponentials to CRAM16/CRAM48.
"""

from unittest.mock import patch, MagicMock

from pytest import approx
import numpy as np
import scipy.sparse as sp
from openmc.deplete.cram import CRAM16, CRAM48
from openmc.deplete.abc import Integrator


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


def test_timed_deplete_clips_negatives():
    """Verify _timed_deplete clips negative atom densities to zero."""
    integrator = MagicMock(spec=Integrator)
    integrator._timed_deplete = Integrator._timed_deplete.__get__(integrator)
    integrator._solver = MagicMock()
    integrator.chain = MagicMock()
    integrator.transfer_rates = None
    integrator.external_source_rates = None

    fake_results = [
        np.array([1.0e10, -1.0e-35, 0.0, -1.0e-47, 5.0e8]),
        np.array([-1.0e-20, 3.0e12, -1.0e-40]),
    ]
    with patch('openmc.deplete.abc.deplete', return_value=fake_results):
        _, results = integrator._timed_deplete(None, None, 1.0)

    for r in results:
        assert np.all(r >= 0.0), f"Negative values found: {r[r < 0]}"
    # Positive values unchanged
    assert results[0][0] == 1.0e10
    assert results[0][4] == 5.0e8
    assert results[1][1] == 3.0e12
    # Negatives clipped to zero
    assert results[0][1] == 0.0
    assert results[0][3] == 0.0
    assert results[1][0] == 0.0
    assert results[1][2] == 0.0
