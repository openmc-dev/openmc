"""Regression tests for openmc.deplete restart capability.

These tests run in two steps, a first run then a restart run, a simple test
problem described in dummy_geometry.py.
"""

import pytest

import openmc.deplete

from tests import dummy_operator


def test_restart_predictor_cecm(run_in_tmpdir):
    """Test to ensure that schemes with different stages are not compatible"""

    op = dummy_operator.DummyOperator()
    output_dir = "test_restart_predictor_cecm"
    op.output_dir = output_dir

    # Perform simulation using the predictor algorithm
    dt = [0.75]
    power = 1.0
    openmc.deplete.PredictorIntegrator(op, dt, power).integrate()

    # Load the files
    prev_res = openmc.deplete.ResultsList.from_hdf5(
        op.output_dir / "depletion_results.h5")

    # Re-create depletion operator and load previous results
    op = dummy_operator.DummyOperator(prev_res)
    op.output_dir = output_dir

    # check ValueError is raised, indicating previous and current stages
    with pytest.raises(ValueError, match="incompatible.* 1.*2"):
        openmc.deplete.CECMIntegrator(op, dt, power)


def test_restart_cecm_predictor(run_in_tmpdir):
    """Integral regression test of integrator algorithm using CE/CM for the
    first run then predictor for the restart run."""

    op = dummy_operator.DummyOperator()
    output_dir = "test_restart_cecm_predictor"
    op.output_dir = output_dir

    # Perform simulation using the MCNPX/MCNP6 algorithm
    dt = [0.75]
    power = 1.0
    cecm = openmc.deplete.CECMIntegrator(op, dt, power)
    cecm.integrate()

    # Load the files
    prev_res = openmc.deplete.ResultsList.from_hdf5(
        op.output_dir / "depletion_results.h5")

    # Re-create depletion operator and load previous results
    op = dummy_operator.DummyOperator(prev_res)
    op.output_dir = output_dir

    # check ValueError is raised, indicating previous and current stages
    with pytest.raises(ValueError, match="incompatible.* 2.*1"):
        openmc.deplete.PredictorIntegrator(op, dt, power)


@pytest.mark.parametrize("scheme", dummy_operator.SCHEMES)
def test_restart(run_in_tmpdir, scheme):
    # set up the problem

    bundle = dummy_operator.SCHEMES[scheme]

    operator = dummy_operator.DummyOperator()

    # take first step
    bundle.solver(operator, [0.75], 1.0).integrate()

    # restart
    prev_res = openmc.deplete.ResultsList.from_hdf5(
        operator.output_dir / "depletion_results.h5")
    operator = dummy_operator.DummyOperator(prev_res)

    # take second step
    bundle.solver(operator, [0.75], 1.0).integrate()

    # compare results

    results = openmc.deplete.ResultsList.from_hdf5(
        operator.output_dir / "depletion_results.h5")

    _t, y1 = results.get_atoms("1", "1")
    _t, y2 = results.get_atoms("1", "2")

    assert y1 == pytest.approx(bundle.atoms_1)
    assert y2 == pytest.approx(bundle.atoms_2)
