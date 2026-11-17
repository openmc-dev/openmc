"""Independent validation; noise-limited bins are explicitly unreliable."""

import numpy as np
from .pce import matrix


def validate(model, unit, observed, errors, *, min_r2=.9, max_noise_mse=2., min_signal_to_noise=3.) -> dict:
    """Calculate per-bin fit quality and conservative sensitivity eligibility."""
    if not np.all(np.isfinite([min_r2,max_noise_mse,min_signal_to_noise])) or min_r2>1 or max_noise_mse<0 or min_signal_to_noise<0:
        raise ValueError("Invalid PCE validation thresholds")
    y = np.asarray(observed, float)
    prediction = model.predict(unit)
    residual = prediction - y
    rmse = np.sqrt(np.mean(residual**2, axis=0))
    spread = np.mean((y - y.mean(axis=0)) ** 2, axis=0)
    r2 = np.divide(
        rmse**2, spread, out=np.full_like(spread, np.inf), where=spread > 0
    )
    r2 = 1 - r2
    a = matrix(unit, model.powers)
    predvar = np.maximum(0, np.einsum("nt,btv,nv->nb", a, model.covariance, a))
    denominator = np.sqrt(np.asarray(errors) ** 2 + predvar)
    standardized = np.divide(
        residual,
        denominator,
        out=np.full_like(residual, np.nan),
        where=denominator > 0,
    )
    noise = np.sqrt(np.mean(np.asarray(errors) ** 2, axis=0))
    coefficient_noise = np.sqrt(
        np.maximum(
            0, np.array([np.trace(c[1:, 1:]) for c in model.covariance])
        )
    )
    signal = np.sqrt(model.variance)
    # Strong raw validation OR noise-consistent residuals, but always require
    # that the inferred parameter spread is resolved above coefficient noise.
    # A closed production channel can have zero observations, uncertainty and
    # surrogate variance. Its standardized residuals are undefined (0/0),
    # so keep an explicit unavailable diagnostic without reducing an empty slice.
    valid_standardized = np.isfinite(standardized)
    noise_samples = valid_standardized.sum(axis=0)
    noise_mse = np.divide(
        np.sum(np.where(valid_standardized, standardized**2, 0.0), axis=0),
        noise_samples,
        out=np.full(prediction.shape[1], np.nan),
        where=noise_samples > 0,
    )
    accepted = (
        ((r2 >= min_r2) | (noise_mse <= max_noise_mse))
        & (signal > min_signal_to_noise * coefficient_noise)
        & np.isfinite(r2)
    )
    accepted[model.diagnostics.get("unavailable_bins", [])] = False
    return dict(
        r2=r2,
        rmse=rmse,
        normalized_rmse=rmse / np.maximum(abs(y.mean(axis=0)), 1e-300),
        maximum_absolute_error=np.max(abs(residual), axis=0),
        noise_mse=noise_mse,
        noise_diagnostic_samples=noise_samples,
        constant_observed_response=spread == 0,
        transport_se_rms=noise,
        signal_to_coefficient_noise=signal
        / np.maximum(coefficient_noise, 1e-300),
        accepted=accepted,
        prediction=prediction,
        residual=residual,
    )
