"""Regression tests for https://github.com/openmc-dev/openmc/issues/3973

`continue_timesteps=True` restarts chained after `integrate(final_step=False)`
were silently zeroing out k/reaction-rate data for every step after the very
first one. Root cause: `Integrator._get_bos_data()` decided whether to reuse
`prev_res[-1]` as beginning-of-step data purely from `step_index == 0`, which
is also true on every restart -- even when `prev_res[-1]` is the zero-filled
placeholder written when a previous `integrate()` call used
`final_step=False`.

These tests exercise the fix at the two points where it lives:

1. `StepResult` now carries an explicit `evaluated` flag through save/load
   (including a default of `True` for pre-existing results files that
   predate the flag).
2. `Integrator._get_bos_data()` now consults that flag -- not just
   `step_index` -- before deciding to reuse restart data.
"""

from unittest.mock import MagicMock

import h5py
import numpy as np
import pytest
from uncertainties import ufloat

from openmc.deplete.abc import Integrator, OperatorResult
from openmc.deplete.stepresult import StepResult


# ---------------------------------------------------------------------------
# StepResult.evaluated: HDF5 round-trip
# ---------------------------------------------------------------------------

class _StubOperator:
    """Just enough of TransportOperator's interface for StepResult.save()"""

    def get_results_info(self):
        volume = {"1": 1.0}
        nuc_list = ["U235"]
        burn_list = ["1"]
        full_burn_list = ["1"]
        name_list = ["fuel"]
        return volume, nuc_list, burn_list, full_burn_list, name_list


def test_stepresult_evaluated_roundtrip(tmp_path):
    """`evaluated` should survive a save -> HDF5 -> load round trip."""
    path = tmp_path / "depletion_results.h5"
    op = _StubOperator()
    x = [np.array([1.0e20])]

    # Step 0: a real, evaluated transport result
    real_result = OperatorResult(k=ufloat(1.0, 0.001), rates=None)
    StepResult.save(op, x, real_result, [0.0, 100.0], 174.0, 0,
                     proc_time=1.0, path=str(path), evaluated=True)

    # Step 1: the zero-filled placeholder written by final_step=False
    placeholder_result = OperatorResult(k=ufloat(0.0, 0.0), rates=None)
    StepResult.save(op, x, placeholder_result, [100.0, 100.0], 174.0, 1,
                     proc_time=1.0, path=str(path), evaluated=False)

    with h5py.File(str(path), "r") as fh:
        step0 = StepResult.from_hdf5(fh, 0)
        step1 = StepResult.from_hdf5(fh, 1)

    assert step0.evaluated is True
    assert step1.evaluated is False


def test_stepresult_evaluated_defaults_true_for_legacy_files(tmp_path):
    """Results files written before this fix (no 'evaluated' dataset)
    should be treated as evaluated=True, preserving old behavior for
    existing files on disk."""
    path = tmp_path / "depletion_results.h5"
    op = _StubOperator()
    x = [np.array([1.0e20])]

    result = OperatorResult(k=ufloat(1.0, 0.001), rates=None)
    StepResult.save(op, x, result, [0.0, 100.0], 174.0, 0, proc_time=1.0,
                     path=str(path))

    # Simulate a results file from before this fix by removing the
    # 'evaluated' dataset that StepResult.save() would normally have
    # written.
    with h5py.File(str(path), "a") as fh:
        del fh["evaluated"]

    with h5py.File(str(path), "r") as fh:
        step0 = StepResult.from_hdf5(fh, 0)

    assert step0.evaluated is True


# ---------------------------------------------------------------------------
# Integrator._get_bos_data(): must not reuse an unevaluated restart step
# ---------------------------------------------------------------------------

class _FakeStepResult:
    """Stand-in for a StepResult; only `.evaluated` is inspected by
    `_get_bos_data`."""

    def __init__(self, evaluated):
        self.evaluated = evaluated


class _FakeOperator:
    def __init__(self, prev_res=None):
        self.prev_res = prev_res


class _MinimalIntegrator(Integrator):
    """Concrete Integrator subclass with the bare minimum to satisfy the
    ABC; `__init__` is never called in these tests (see `_make_integrator`
    below), so its signature/body doesn't matter here."""

    _num_stages = 1

    def __call__(self, n, rates, dt, source_rate, i):
        raise NotImplementedError


def _make_integrator(prev_res):
    """Build an Integrator instance without running __init__ (which would
    require a real operator, chain, timesteps, etc.), and stub out the two
    branches that `_get_bos_data` chooses between so we can observe which
    one gets called."""
    integ = object.__new__(_MinimalIntegrator)
    integ.operator = _FakeOperator(prev_res=prev_res)
    integ._keff_search_control = None
    integ._i_res = 0
    integ._get_bos_data_from_operator = MagicMock(
        return_value=("bos_conc", "from_operator"))
    integ._get_bos_data_from_restart = MagicMock(
        return_value=("bos_conc", "from_restart"))
    return integ


def test_get_bos_data_uses_operator_when_no_prev_results():
    integ = _make_integrator(prev_res=None)
    integ._get_bos_data(step_index=0, source_rate=174.0, bos_conc=None)
    integ._get_bos_data_from_operator.assert_called_once()
    integ._get_bos_data_from_restart.assert_not_called()


def test_get_bos_data_uses_operator_when_step_index_nonzero():
    integ = _make_integrator(prev_res=[_FakeStepResult(evaluated=True)])
    integ._get_bos_data(step_index=1, source_rate=174.0, bos_conc=None)
    integ._get_bos_data_from_operator.assert_called_once()
    integ._get_bos_data_from_restart.assert_not_called()


def test_get_bos_data_reuses_restart_when_evaluated():
    """Normal restart case: continuing from a run that ended with
    final_step=True (or the natural end of the loop) should still reuse
    the real, previously-computed data -- this must keep working."""
    integ = _make_integrator(prev_res=[_FakeStepResult(evaluated=True)])
    integ._get_bos_data(step_index=0, source_rate=174.0, bos_conc=None)
    integ._get_bos_data_from_restart.assert_called_once()
    integ._get_bos_data_from_operator.assert_not_called()


def test_get_bos_data_ignores_unevaluated_restart_step():
    """Regression test for #3973: when the last StepResult in prev_res is
    the zero-filled final_step=False placeholder, _get_bos_data must NOT
    treat it as valid restart data -- it must fall back to a real
    transport evaluation instead."""
    integ = _make_integrator(prev_res=[_FakeStepResult(evaluated=False)])
    integ._get_bos_data(step_index=0, source_rate=174.0, bos_conc=None)
    integ._get_bos_data_from_operator.assert_called_once()
    integ._get_bos_data_from_restart.assert_not_called()


# ---------------------------------------------------------------------------
# Integrator.integrate(): must not write a statepoint for a skipped
# (final_step=False) trailing evaluation
# ---------------------------------------------------------------------------

class _FakeOperatorForIntegrate:
    """Stub TransportOperator sufficient to drive Integrator.integrate()
    through the trailing 'final simulation' block without touching real
    OpenMC/transport machinery."""

    def __init__(self):
        self.prev_res = None
        self.output_dir = "."
        self.chain = None
        self.write_bos_data = MagicMock()
        self._call_source_rates = []

    def initial_condition(self):
        return [np.array([1.0e20])]

    def get_results_info(self):
        volume = {"1": 1.0}
        nuc_list = ["U235"]
        burn_list = ["1"]
        full_burn_list = ["1"]
        name_list = ["fuel"]
        return volume, nuc_list, burn_list, full_burn_list, name_list

    def __call__(self, vec, source_rate):
        # Mirrors CoupledOperator.__call__: a zero source rate means the
        # operator was reset() but never actually run().
        self._call_source_rates.append(source_rate)
        if source_rate == 0.0:
            return OperatorResult(k=ufloat(0.0, 0.0), rates=None)
        return OperatorResult(k=ufloat(1.15, 0.003), rates=None)

    def finalize(self):
        pass


class _OneStepIntegrator(Integrator):
    """Minimal concrete Integrator: one Bateman step that just returns the
    input concentrations unchanged (no real depletion math needed to
    exercise integrate()'s control flow)."""

    _num_stages = 1

    def __call__(self, n, rates, dt, source_rate, i):
        return 0.0, n


def _run_integrate(tmp_path, final_step):
    op = _FakeOperatorForIntegrate()
    integ = _OneStepIntegrator(
        op, timesteps=[1.0], power=174.0, timestep_units="d")
    integ.integrate(
        final_step=final_step, output=False,
        path=str(tmp_path / "depletion_results.h5"))
    return op


def test_integrate_skips_write_bos_data_when_final_step_false(tmp_path):
    """Regression test: with final_step=False, the trailing 'final
    simulation' evaluation never actually runs transport (source_rate=0
    short-circuits it), so write_bos_data() -- which triggers a
    statepoint write -- must not be called for that step. Previously it
    was called unconditionally, producing statepoints with leftover
    internal state (observed as a finite k-effective paired with an
    infinite standard deviation)."""
    op = _run_integrate(tmp_path, final_step=False)
    # One real call for step 0 (source_rate=174.0), one skipped call for
    # the trailing placeholder (source_rate=0.0)
    assert op._call_source_rates == [174.0, 0.0]
    # write_bos_data should only have been invoked for the real step 0
    # evaluation (from _get_bos_data_from_operator), not for the skipped
    # trailing placeholder.
    op.write_bos_data.assert_called_once_with(0)


def test_integrate_calls_write_bos_data_when_final_step_true(tmp_path):
    """With final_step=True, the trailing evaluation is real, so its
    statepoint should still be written as before."""
    op = _run_integrate(tmp_path, final_step=True)
    assert op._call_source_rates == [174.0, 174.0]
    assert op.write_bos_data.call_count == 2
    op.write_bos_data.assert_any_call(0)
    op.write_bos_data.assert_any_call(1)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))
