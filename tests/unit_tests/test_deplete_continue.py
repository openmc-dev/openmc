"""Unit tests for openmc.deplete continue run capability.

These tests run in two steps: first a normal run and then a continue run based on the prev_results
"""

import pytest

import openmc.deplete

from tests import dummy_operator

# test that the continue timesteps works when the second integrate call contains all previous timesteps
@pytest.mark.parametrize("scheme", dummy_operator.SCHEMES)
def test_continue(run_in_tmpdir, scheme):
    # set up the problem

    bundle = dummy_operator.SCHEMES[scheme]

    operator = dummy_operator.DummyOperator()

    # take first step
    bundle.solver(operator, [0.75], 1.0).integrate()

    # restart
    prev_res = openmc.deplete.Results(
        operator.output_dir / "depletion_results.h5")
    operator = dummy_operator.DummyOperator(prev_res)

    # if continue run happens, test passes
    bundle.solver(operator, [0.75, 0.75], [1.0, 1.0], continue_timesteps = True).integrate()

@pytest.mark.parametrize("scheme", dummy_operator.SCHEMES)
def test_mismatched_initial_steps(run_in_tmpdir, scheme):
    """Test to ensure that a continue run with different initial steps is properly caught"""

    # set up the problem

    bundle = dummy_operator.SCHEMES[scheme]

    operator = dummy_operator.DummyOperator()

    # take first step
    bundle.solver(operator, [0.75, 0.75], [1.0,1.0]).integrate()

    # restart
    prev_res = openmc.deplete.Results(
        operator.output_dir / "depletion_results.h5")
    operator = dummy_operator.DummyOperator(prev_res)

    # continue run with different previous step should cause a ValueError
    # note the first step matches but the second does not while the third is a new step
    with pytest.raises(ValueError):
        bundle.solver(operator, [0.75, 0.5, 0.75], [1.0, 2.0, 1.0], continue_timesteps = True).integrate()