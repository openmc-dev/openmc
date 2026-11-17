"""Exact variance decomposition of orthonormal PCE coefficients."""

import numpy as np


def analyse(model, groups=None, accepted=None) -> dict:
    """Return main, total, pure pair and closed/group-total indices.

    Group first order includes terms supported entirely within the group;
    group total includes every term touching it. Individual totals are NOT
    summed. Zero variance and rejected bins yield NaN, never confident zeros.
    """
    support = model.powers > 0
    power = model.coefficients**2
    power[0] = 0
    variance = power.sum(axis=0)
    valid = variance > 0
    if accepted is not None:
        valid &= accepted

    def fraction(mask):
        return np.divide(
            power[mask].sum(axis=0),
            variance,
            out=np.full_like(variance, np.nan),
            where=valid,
        )

    first = np.array(
        [
            fraction(support[:, i] & (support.sum(axis=1) == 1))
            for i in range(support.shape[1])
        ]
    )
    total = np.array(
        [fraction(support[:, i]) for i in range(support.shape[1])]
    )
    second = {
        f"{i},{j}": fraction(
            support[:, i] & support[:, j] & (support.sum(axis=1) == 2)
        )
        for i in range(support.shape[1])
        for j in range(i + 1, support.shape[1])
    }
    diagonal = np.diagonal(model.covariance, axis1=1, axis2=2)

    def signal_to_noise(mask):
        noise = np.maximum(0, diagonal[:, mask].sum(axis=1))
        return np.sqrt(power[mask].sum(axis=0) / np.maximum(noise, 1e-300))

    component_snr = np.array(
        [signal_to_noise(support[:, i]) for i in range(support.shape[1])]
    )
    grouped = {}
    for name, members in (groups or {}).items():
        touching = support[:, members].any(axis=1)
        outside = np.delete(support, members, axis=1).any(axis=1)
        grouped[name] = dict(
            first=fraction(touching & ~outside),
            total=fraction(touching),
            signal_to_coefficient_noise=signal_to_noise(touching),
            resolved=(signal_to_noise(touching) > 3) & valid,
        )
    return dict(
        first=first,
        total=total,
        second=second,
        groups=grouped,
        variance=variance,
        total_component_signal_to_noise=component_snr,
        total_component_resolved=(component_snr > 3) & valid[None, :],
    )


def stability(model, groups=None, draws=200, seed=97) -> dict:
    """Refit-noise propagation intervals; conditional on basis/truncation."""
    from .pce import PCE

    rng = np.random.default_rng(seed)
    samples = []
    for _ in range(draws):
        c = np.stack(
            [
                rng.multivariate_normal(
                    model.coefficients[:, j],
                    model.covariance[j],
                    check_valid="ignore",
                )
                for j in range(model.coefficients.shape[1])
            ],
            axis=1,
        )
        samples.append(
            analyse(PCE(model.powers, c, model.covariance, {}), groups)[
                "total"
            ]
        )
    return dict(
        total_q05=np.quantile(samples, 0.05, axis=0),
        total_q95=np.quantile(samples, 0.95, axis=0),
    )
