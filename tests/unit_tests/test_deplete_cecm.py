"""Regression tests for openmc.deplete.integrator.cecm algorithm.

These tests integrate a simple test problem described in dummy_geometry.py.
"""

from pytest import approx
import openmc.deplete
from openmc.deplete import results
from openmc.deplete import utilities

from tests import dummy_geometry


def test_cecm(run_in_tmpdir):
    """Integral regression test of integrator algorithm using CE/CM."""

    settings = openmc.deplete.Settings()
    settings.dt_vec = [0.75, 0.75]
    settings.output_dir = "test_integrator_regression"

    op = dummy_geometry.DummyGeometry(settings)

    # Perform simulation using the MCNPX/MCNP6 algorithm
    openmc.deplete.cecm(op, print_out=False)

    # Load the files
    res = results.read_results(settings.output_dir + "/results.h5")

    _, y1 = utilities.evaluate_single_nuclide(res, "1", "1")
    _, y2 = utilities.evaluate_single_nuclide(res, "1", "2")

    # Mathematica solution
    s1 = [1.86872629872102, 1.395525772416039]
    s2 = [2.18097439443550, 2.69429754646747]

    assert y1[1] == approx(s1[0])
    assert y2[1] == approx(s1[1])

    assert y1[2] == approx(s2[0])
    assert y2[2] == approx(s2[1])
