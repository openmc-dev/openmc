"""Genuine Morris grid trajectories and noise-aware trajectory statistics.

Coordinates are standardized to [0,1]. Effects are response change divided by
signed unit-coordinate step; physical coefficient step is 2*limit*delta.
"""

from dataclasses import dataclass
import numpy as np
from scipy.stats import spearmanr


@dataclass
class MorrisDesign:
    points: np.ndarray
    changed: np.ndarray
    steps: np.ndarray
    levels: int
    seed: int


def trajectories(
    dimension: int, count: int, *, levels=6, seed=2026
) -> MorrisDesign:
    """Generate randomized OAT trajectories using delta=p/[2(p−1)].

    Each point lies on the p-level grid. Each coordinate is changed once per
    trajectory, in a randomized order, with a valid random step orientation.
    Increasing count preserves the prefix for exact resumability.
    """
    if dimension < 1 or count < 1 or levels < 4 or levels % 2:
        raise ValueError(
            "Positive dimensions/count and even levels>=4 required"
        )
    rng = np.random.default_rng(seed)
    jump = levels // 2
    points = []
    changed = []
    steps = []
    for _ in range(count):
        low = rng.integers(0, levels - jump, size=dimension)
        direction = rng.choice([-1, 1], size=dimension)
        start = low + np.where(direction < 0, jump, 0)
        order = rng.permutation(dimension)
        path = [start.copy()]
        signed = []
        for j in order:
            start = start.copy()
            start[j] += direction[j] * jump
            path.append(start)
            signed.append(direction[j] * jump / (levels - 1))
        points.append(np.array(path) / (levels - 1))
        changed.append(order)
        steps.append(signed)
    return MorrisDesign(
        np.array(points), np.array(changed), np.array(steps), levels, seed
    )


def effects(design: MorrisDesign, means: np.ndarray, errors: np.ndarray):
    """Return EE and independent-endpoint EE SE, shape (r,d,bins).

    Adjacent effects share observations; resampling entire trajectories retains
    this dependence. Standard errors assume independent transport seeds at
    trajectory vertices; same-seed experiments require measured covariance.
    """
    r, n, d = design.points.shape
    means = np.asarray(means)
    errors = np.asarray(errors)
    if means.shape[:2] != (r, n) or means.shape != errors.shape:
        raise ValueError("Response arrays do not match Morris design")
    ee = np.empty((r, d, means.shape[2]))
    se = np.empty_like(ee)
    for t in range(r):
        for k, j in enumerate(design.changed[t]):
            step = design.steps[t, k]
            ee[t, j] = (means[t, k + 1] - means[t, k]) / step
            se[t, j] = np.sqrt(
                errors[t, k + 1] ** 2 + errors[t, k] ** 2
            ) / abs(step)
    return ee, se


def statistics(ee, se, *, bootstraps=500, seed=19, confidence_level=.95, noise_multiplier=2.):
    """Bootstrap complete trajectories, report noise and convergence separately."""
    if not 0 < confidence_level < 1 or bootstraps < 1 or not np.isfinite(noise_multiplier) or noise_multiplier < 0:
        raise ValueError("Invalid Morris bootstrap or noise settings")
    ee = np.asarray(ee)
    se = np.asarray(se)
    r = ee.shape[0]
    if r < 2:
        raise ValueError("At least two trajectories required for variability")
    mu = ee.mean(axis=0)
    star = abs(ee).mean(axis=0)
    sigma = ee.std(axis=0, ddof=1)
    noise = np.sqrt(np.mean(se**2, axis=0))
    # Signed mean SE gives a conservative resolution scale; mu* is positive
    # under transport noise, so it is never treated as automatically resolved.
    mean_noise = np.sqrt(np.sum(se**2, axis=0)) / r
    snr = star / np.maximum(noise, 1e-300)
    rng = np.random.default_rng(seed)
    resamples = rng.integers(0, r, size=(bootstraps, r))
    low = np.empty_like(star)
    high = np.empty_like(star)
    mu_low = np.empty_like(mu)
    mu_high = np.empty_like(mu)
    sigma_low = np.empty_like(sigma)
    sigma_high = np.empty_like(sigma)
    # Bound spatial-bootstrap memory; never allocate (bootstrap,r,d,all_bins).
    for start in range(0, ee.shape[2], 16):
        stop = start + 16
        part = ee[:, :, start:stop]
        for metric, left_ci, right_ci in (
            ("star", low, high),
            ("mu", mu_low, mu_high),
            ("sigma", sigma_low, sigma_high),
        ):
            draws = np.array(
                [
                    (
                        abs(part[idx]).mean(axis=0)
                        if metric == "star"
                        else (
                            part[idx].mean(axis=0)
                            if metric == "mu"
                            else part[idx].std(axis=0, ddof=1)
                        )
                    )
                    for idx in resamples
                ]
            )
            bounds = np.quantile(draws, [(1-confidence_level)/2, (1+confidence_level)/2], axis=0)
            left_ci[:, start:stop] = bounds[0]
            right_ci[:, start:stop] = bounds[1]
    correction = np.sqrt(
        np.maximum(0, np.mean(ee**2, axis=0) - np.mean(se**2, axis=0))
    )
    correction[~np.isfinite(correction)] = np.nan
    half = r // 2
    left = abs(ee[:half]).mean(axis=0)
    right = abs(ee[half:]).mean(axis=0)
    rank_correlation = []
    overlap = []
    for b in range(ee.shape[2]):
        k = min(10, ee.shape[1])
        rank_correlation.append(
            float(spearmanr(left[:, b], right[:, b]).statistic)
            if ee.shape[1] > 1
            and np.ptp(left[:, b]) > 0
            and np.ptp(right[:, b]) > 0
            else (1.0 if ee.shape[1] == 1 else float("nan"))
        )
        overlap.append(
            len(
                set(np.argsort(-left[:, b])[:k])
                & set(np.argsort(-right[:, b])[:k])
            )
            / k
        )
    return dict(
        mu=mu,
        mu_star=star,
        sigma=sigma,
        mu_ci_low=mu_low,
        mu_ci_high=mu_high,
        sigma_ci_low=sigma_low,
        sigma_ci_high=sigma_high,
        mu_star_ci_low=low,
        mu_star_ci_high=high,
        noise_rms=noise,
        mean_noise_se=mean_noise,
        signal_to_noise=snr,
        noise_corrected_rms=correction,
        resolved=(low > noise_multiplier * noise) & (snr > noise_multiplier),
        sigma_excess=np.sqrt(np.maximum(0, sigma**2 - noise**2)),
        ranking_spearman_halves=np.array(rank_correlation),
        top10_overlap=np.array(overlap),
        trajectories=r,
        bootstrap_unit="whole trajectory",
    )


def convergence(ee, se, minimum=20, *, rank_correlation=.9, top_overlap=.8,
                max_relative_ci_width=.5, bootstraps=200, seed=19,
                confidence_level=.95, noise_multiplier=2.):
    """Assess ranks as r grows; stopping requires >=20 trajectories and precision.

    Stable top-ten overlap>=.8, rank rho>=.9, and median resolved top-ten
    relative CI width<=.5 on the total response. This is a declared operational
    criterion, not a mathematical guarantee of global ranking convergence.
    """
    rows = []
    for count in sorted({2, 4, 8, 12, 20, 40, 60, len(ee)}):
        if count > len(ee):
            continue
        s = statistics(ee[:count], se[:count], bootstraps=bootstraps, seed=seed, confidence_level=confidence_level, noise_multiplier=noise_multiplier)
        top = np.argsort(-s["mu_star"][:, 0])[: min(10, ee.shape[1])]
        width = np.median(
            (s["mu_star_ci_high"][top, 0] - s["mu_star_ci_low"][top, 0])
            / np.maximum(s["mu_star"][top, 0], 1e-300)
        )
        stopped = bool(
            count >= minimum
            and s["ranking_spearman_halves"][0] >= rank_correlation
            and s["top10_overlap"][0] >= top_overlap
            and width <= max_relative_ci_width
            and s["resolved"][top, 0].all()
        )
        rows.append(
            dict(
                trajectories=count,
                rho=float(s["ranking_spearman_halves"][0]),
                top10_overlap=float(s["top10_overlap"][0]),
                median_relative_ci_width=float(width),
                stopping_criterion_met=stopped,
            )
        )
    return rows
