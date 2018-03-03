"""Regression tests for openmc.deplete.integrator.predictor algorithm.

These tests integrate a simple test problem described in dummy_geometry.py.
"""

from pytest import approx
import openmc.deplete

from tests import dummy_operator


def test_predictor(run_in_tmpdir):
    """Integral regression test of integrator algorithm using predictor/corrector"""

    op = dummy_operator.DummyOperator()
    op.output_dir = "test_integrator_regression"

    # Perform simulation using the predictor algorithm
    dt = [0.75, 0.75]
    power = 1.0
    openmc.deplete.predictor(op, dt, power, print_out=False)

    # Load the files
    res = openmc.deplete.ResultsList(op.output_dir / "depletion_results.h5")

    _, y1 = res.get_atoms("1", "1")
    _, y2 = res.get_atoms("1", "2")

    # Mathematica solution
    s1 = [2.46847546272295, 0.986431226850467]
    s2 = [4.11525874568034, -0.0581692232513460]

    assert y1[1] == approx(s1[0])
    assert y2[1] == approx(s1[1])

    assert y1[2] == approx(s2[0])
    assert y2[2] == approx(s2[1])
