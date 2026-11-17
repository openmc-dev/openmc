"""Deterministic validation of the example's UQ math, separate from transport."""

import numpy as np
from numpy.polynomial.legendre import leggauss
import pytest

from sandbox.examples.fusion_li6_uq_scan import (
    multi_indices,
    pce_matrix,
    fit_pce,
    predict,
    field_statistics,
)


def test_legendre_product_orthonormality():
    x, weights = leggauss(3)
    points = np.array([(a, b) for a in (x + 1) / 2 for b in (x + 1) / 2])
    joint_weights = np.array([a * b / 4 for a in weights for b in weights])
    indices = multi_indices(2, 2)
    matrix = pce_matrix(points, indices)
    np.testing.assert_allclose(
        matrix.T @ (joint_weights[:, None] * matrix),
        np.eye(len(indices)),
        atol=2.0e-15,
    )


def test_exact_log_polynomial_recovered():
    rng = np.random.default_rng(19)
    training = rng.random((64, 2))
    validation = rng.random((32, 2))

    def response(u):
        return np.exp(
            0.5 + 0.2 * u[:, 0] - 0.3 * u[:, 1] ** 2 + 0.1 * u[:, 0] * u[:, 1]
        )

    rates = response(training)
    indices, coefficients, covariance = fit_pce(
        training, rates, 0.001 * rates, 2
    )
    np.testing.assert_allclose(
        predict(validation, indices, coefficients),
        response(validation),
        rtol=2.0e-14,
    )
    assert np.linalg.eigvalsh(covariance).min() > 0


def test_rank_deficiency_rejected():
    with pytest.raises(RuntimeError, match="Rank-deficient"):
        fit_pce(np.zeros((64, 2)), np.ones(64), np.full(64, 0.01), 2)


def test_transport_noise_is_not_reported_as_parameter_variance():
    # Observed sample variance=5; average transport variance=1, so corrected=4.
    means = np.array([[1.0, 1.0], [2.0, 2.0], [4.0, 4.0], [6.0, 6.0]])
    errors = np.ones_like(means)
    errors[0, 1] = np.inf
    average, corrected, noise, variance = field_statistics(means, errors)
    np.testing.assert_allclose(average, [3.25, 3.25])
    np.testing.assert_allclose(corrected[0] ** 2, variance[0] - 1.0)
    assert noise[0] == 1.0
    assert np.isnan(corrected[1])


def test_pce_only_propagation_matches_log_uniform_distribution(
    tmp_path, monkeypatch
):
    import json
    import openmc.lib
    from sandbox.examples.propagate_li6_pce import propagate

    # A known single-Legendre-term log-rate model has an exact log-uniform
    # response. It supplies an independent reference for propagated moments,
    # quantiles and coefficient-noise propagation without native transport.
    directory = tmp_path / "model"
    directory.mkdir()
    np.savez(
        directory / "pce_model.npz",
        multi_indices=np.array([[0, 0, 0, 0, 0], [1, 0, 0, 0, 0]]),
        coefficients=np.array([np.log(0.16), 0.05]),
        coefficient_covariance_from_tally_noise=np.diag([1.0e-8, 0.0]),
    )
    (directory / "metadata.json").write_text(
        json.dumps(
            dict(
                pce=dict(
                    validation=dict(accepted=True),
                    predicted_distribution=dict(baseline_mean=0.16),
                )
            )
        )
    )

    def no_transport(*args, **kwargs):
        pytest.fail("Saved PCE propagation must not initialize or run OpenMC")

    monkeypatch.setattr(openmc.lib, "init", no_transport)
    monkeypatch.setattr(openmc.lib, "run", no_transport)
    summary = propagate(directory, 10000, 29, tmp_path / "output")
    half_width = 0.05 * np.sqrt(3.0)
    exact_mean = 0.16 * np.sinh(half_width) / half_width
    assert (
        abs(summary["mean"] - exact_mean)
        < 5 * summary["mean_sampling_standard_error"]
    )
    exact_quantiles = 0.16 * np.exp(
        half_width * (2 * np.array([0.05, 0.5, 0.95]) - 1)
    )
    errors = np.array(summary["quantile_sampling_standard_errors"])
    assert np.all(
        abs(np.array(summary["quantiles_5_50_95"]) - exact_quantiles)
        < 5 * errors
    )
    assert summary[
        "mean_fit_standard_error_from_training_tally_noise"
    ] == pytest.approx(summary["mean"] * 1.0e-4)
