"""Morris validity and native-node energy partition invariants."""

import numpy as np
import pytest
from sandbox.analysis_codes.morris import (
    trajectories,
    effects,
    statistics,
    convergence,
)
from sandbox.analysis_codes.registry import Parameter
from sandbox.analysis_codes.cross_sections import perturb


@pytest.mark.parametrize("levels", [6, 8])
def test_grid_trajectory_validity(levels):
    d = 7
    r = 20
    p = trajectories(d, r, levels=levels, seed=12)
    assert p.points.shape == (r, d + 1, d)
    assert np.all((p.points >= 0) & (p.points <= 1))
    np.testing.assert_allclose(
        p.points * (levels - 1), np.round(p.points * (levels - 1)), atol=1e-14
    )
    np.testing.assert_array_equal(
        np.sort(p.changed, axis=1), np.tile(np.arange(d), (r, 1))
    )
    assert np.all(np.count_nonzero(np.diff(p.points, axis=1), axis=2) == 1)
    np.testing.assert_allclose(abs(p.steps), levels / (2 * (levels - 1)))
    np.testing.assert_array_equal(
        trajectories(d, 40, levels=levels, seed=12).points[:r], p.points
    )


def test_analytic_elementary_effects():
    p = trajectories(3, 40, seed=8)
    slopes = np.array([1.0, -2.0, 3.0])
    y = (5 + p.points @ slopes)[..., None]
    errors = np.zeros_like(y)
    ee, se = effects(p, y, errors)
    s = statistics(ee, se, bootstraps=200)
    np.testing.assert_allclose(s["mu"][:, 0], slopes, atol=1e-14)
    np.testing.assert_allclose(s["mu_star"][:, 0], abs(slopes), atol=1e-14)
    np.testing.assert_allclose(s["sigma"], 0, atol=1e-14)
    assert convergence(ee, se)[-1]["stopping_criterion_met"]


def test_morris_independent_noise():
    p = trajectories(1, 4, seed=1)
    means = np.zeros((4, 2, 1))
    errors = np.ones_like(means) * 0.1
    ee, se = effects(p, means, errors)
    np.testing.assert_allclose(se, np.sqrt(2) * 0.1 / abs(p.steps)[:, :, None])
    assert not statistics(ee, se)["resolved"].any()




def test_interacting_nonlinear_elementary_effects():
    plan = trajectories(3, 20, seed=812)
    x = plan.points
    y = (x[:, :, 0] * x[:, :, 1] + x[:, :, 2] ** 2)[..., None]
    ee, se = effects(plan, y, np.zeros_like(y))
    for t in range(len(x)):
        for k, j in enumerate(plan.changed[t]):
            expected = (
                x[t, k, 1]
                if j == 0
                else (
                    x[t, k, 0] if j == 1 else 2 * x[t, k, 2] + plan.steps[t, k]
                )
            )
            assert ee[t, j, 0] == pytest.approx(expected)
    assert statistics(ee, se)["sigma_excess"].max() > 0
