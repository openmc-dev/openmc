"""Unit tests for openmc.deplete continue run capability.

These tests run in two steps: first a normal run and then a continue run based on the prev_results
"""

import pytest

import openmc.deplete

from tests import dummy_operator

@pytest.mark.parametrize("scheme", dummy_operator.SCHEMES)
def test_continue(run_in_tmpdir, scheme):
    """Test to ensure that a properly defined continue run works"""

    # set up the problem

    bundle = dummy_operator.SCHEMES[scheme]

    operator = dummy_operator.DummyOperator()

    # initial depletion
    bundle.solver(operator, [1.0, 2.0], [1.0, 2.0]).integrate()

    # set up continue run
    prev_res = openmc.deplete.Results(operator.output_dir / "depletion_results.h5")
    operator = dummy_operator.DummyOperator(prev_res)

    # if continue run happens, test passes
    bundle.solver(operator, [1.0, 2.0, 3.0, 4.0], [1.0, 2.0, 3.0, 4.0], continue_timesteps=True).integrate()

@pytest.mark.parametrize("scheme", dummy_operator.SCHEMES)
def test_continue_for_null_previous(run_in_tmpdir, scheme):
    """Test to ensure that a continue run works even if there are no previous results"""
    # set up the problem

    bundle = dummy_operator.SCHEMES[scheme]

    operator = dummy_operator.DummyOperator()

    # initial depletion
    bundle.solver(operator, [1.0, 2.0], [1.0, 2.0], continue_timesteps=True).integrate()

@pytest.mark.parametrize("scheme", dummy_operator.SCHEMES)
def test_mismatched_initial_times(run_in_tmpdir, scheme):
    """Test to ensure that a continue run with different initial steps is properly caught"""

    # set up the problem

    bundle = dummy_operator.SCHEMES[scheme]

    operator = dummy_operator.DummyOperator()

    # perform initial steps
    bundle.solver(operator, [0.75, 0.75], [1.0, 1.0]).integrate()

    # restart
    prev_res = openmc.deplete.Results(operator.output_dir / "depletion_results.h5")
    operator = dummy_operator.DummyOperator(prev_res)

    with pytest.raises(
        ValueError,
        match="You are attempting to continue a run in which the previous timesteps "
        "do not have the same initial timesteps as those provided to the "
        "Integrator. Please make sure you are using the correct timesteps.",
    ):
        bundle.solver(
            operator, [0.75, 0.5, 0.75], [1.0, 1.0, 1.0], continue_timesteps=True
        ).integrate()


@pytest.mark.parametrize("scheme", dummy_operator.SCHEMES)
def test_mismatched_initial_source_rates(run_in_tmpdir, scheme):
    """Test to ensure that a continue run with different initial steps is properly caught"""

    # set up the problem

    bundle = dummy_operator.SCHEMES[scheme]

    operator = dummy_operator.DummyOperator()

    # perform initial steps
    bundle.solver(operator, [0.75, 0.75], [1.0, 1.0]).integrate()

    # restart
    prev_res = openmc.deplete.Results(operator.output_dir / "depletion_results.h5")
    operator = dummy_operator.DummyOperator(prev_res)

    with pytest.raises(
        ValueError,
        match="You are attempting to continue a run in which the previous results "
        "do not have the same initial source rates, powers, or power densities "
        "as those provided to the Integrator. Please make sure you are using "
        "the correct powers, power densities, or source rates and previous results file.",
    ):
        bundle.solver(
            operator, [0.75, 0.75, 0.75], [1.0, 2.0, 1.0], continue_timesteps=True
        ).integrate()
