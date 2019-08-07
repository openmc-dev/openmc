"""Regression tests for openmc.deplete.integrator.cecm algorithm.

These tests integrate a simple test problem described in dummy_geometry.py.
"""

from pytest import approx
from openmc.deplete import CECMIntegrator, ResultsList

from tests import dummy_operator


def test_cecm(run_in_tmpdir):
    """Integral regression test of integrator algorithm using CE/CM."""

    op = dummy_operator.DummyOperator()
    op.output_dir = "test_integrator_regression"

    # Perform simulation using the MCNPX/MCNP6 algorithm
    dt = [0.75, 0.75]
    power = 1.0
    integrator = CECMIntegrator(op, dt, power)
    integrator.integrate()

    # Load the files
    res = ResultsList.from_hdf5(op.output_dir / "depletion_results.h5")

    _, y1 = res.get_atoms("1", "1")
    _, y2 = res.get_atoms("1", "2")

    # Mathematica solution
    s1 = [1.86872629872102, 1.395525772416039]
    s2 = [2.18097439443550, 2.69429754646747]

    assert y1[1] == approx(s1[0])
    assert y2[1] == approx(s1[1])

    assert y1[2] == approx(s2[0])
    assert y2[2] == approx(s2[1])

    # Test structure of depletion time dataset

    dep_time = res.get_depletion_time()
    assert dep_time.shape == (len(dt), )
    assert all(dep_time > 0)
