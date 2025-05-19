"""Unit tests for openmc.deplete continue run capability.

These tests run in two steps: first a normal run and then a continue run using the previous results
"""

import pytest
import numpy as np
import openmc.deplete
from pathlib import Path

from tests import dummy_operator


def test_continue(run_in_tmpdir):
    """Test to ensure that a properly defined continue run works"""
    # set up the problem
    bundle = dummy_operator.SCHEMES['predictor']
    operator = dummy_operator.DummyOperator()

    # initial depletion
    bundle.solver(operator, [1.0, 2.0], [1.0, 2.0]).integrate()

    # set up continue run
    prev_res = openmc.deplete.Results(operator.output_dir / "depletion_results.h5")
    operator = dummy_operator.DummyOperator(prev_res)

    # if continue run happens, test passes
    bundle.solver(operator, [1.0, 2.0, 3.0, 4.0], [1.0, 2.0, 3.0, 4.0],
                  continue_timesteps=True).integrate()

    final_res = openmc.deplete.Results(operator.output_dir / "depletion_results.h5")

    assert np.array_equal(np.diff(final_res.get_times(time_units="s")),[1.0, 2.0, 3.0, 4.0])


def test_continue_continue(run_in_tmpdir):
    """Test to ensure that a continue run can be continued"""
    # set up the problem
    bundle = dummy_operator.SCHEMES['predictor']
    operator = dummy_operator.DummyOperator()

    # initial depletion
    bundle.solver(operator, [1.0, 2.0], [1.0, 2.0]).integrate()

    # set up continue run
    prev_res = openmc.deplete.Results(operator.output_dir / "depletion_results.h5")
    operator = dummy_operator.DummyOperator(prev_res)

    # first continue run
    bundle.solver(operator, [1.0, 2.0, 3.0, 4.0], [1.0, 2.0, 3.0, 4.0],
                  continue_timesteps=True).integrate()

    prev_res = openmc.deplete.Results(operator.output_dir / "depletion_results.h5")
    # second continue run
    bundle.solver(operator, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0], [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
                  continue_timesteps=True).integrate()
    
    final_res = openmc.deplete.Results(operator.output_dir / "depletion_results.h5")

    assert np.array_equal(np.diff(final_res.get_times(time_units="s")),[1.0, 2.0, 3.0, 4.0, 5.0, 6.0])

def test_killed_and_continue(run_in_tmpdir):
    """
    Attempt to continue from a simulation that was killed mid state.
    The previous state is provided in the form of a few local files:

    continue_model.xml contains the necessary XML information
    simple_chain.xml contians a simplified version of the CASL chain

    continue_depletion_results.h5 contains the results output by
    an OpenMC (v0.15.2) depletion simulation that was killed in
    the middle of the third step
    """
    base_path = Path(__file__).parents[0]
    model = openmc.Model.from_model_xml(f"{base_path}/kill_continue/continue_model.xml")
    power = 35000 # W

    time_steps = [1.0,2.0,3.0,4.0] # days
    prev_results = openmc.deplete.Results(f"{base_path}/kill_continue/continue_depletion_results.h5")
    operator = openmc.deplete.CoupledOperator(
        model, prev_results=prev_results, chain_file=f"{base_path}/kill_continue/simple_chain.xml"
    )
    integrator = openmc.deplete.CECMIntegrator(
        operator,
        time_steps,
        power=power,
        timestep_units="d",
        continue_timesteps=True,
    )
    integrator.integrate(path=f"{base_path}/kill_continue/continue_depletion_results.h5")
    final_res = openmc.deplete.Results(f"{base_path}/kill_continue/continue_depletion_results.h5")

    assert np.array_equal(np.diff(final_res.get_times(time_units="d")),[1.0, 2.0, 3.0, 4.0])

def test_mismatched_initial_times(run_in_tmpdir):
    """Test to ensure that a continue run with different initial steps is properly caught"""
    # set up the problem
    bundle = dummy_operator.SCHEMES['predictor']
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


def test_mismatched_initial_source_rates(run_in_tmpdir):
    """Test to ensure that a continue run with different initial steps is properly caught"""
    # set up the problem
    bundle = dummy_operator.SCHEMES['predictor']
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
