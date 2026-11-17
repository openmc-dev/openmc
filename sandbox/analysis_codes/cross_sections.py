"""Threshold-preserving microscopic cross-section candidates."""

import numpy as np


def perturb(
    values: np.ndarray, grid: np.ndarray, parameters, coefficients
) -> np.ndarray:
    """Apply disjoint linear or log-fractional incident-energy parameters."""
    values = np.asarray(values, float)
    if not np.all(np.isfinite(values)) or np.any(values < 0):
        raise ValueError("Invalid baseline cross sections")
    change = np.zeros_like(values)
    coverage = np.zeros(len(values), int)
    modes = {p.xs_mode for p in parameters}
    if len(modes) != 1:
        raise ValueError("Mixed XS parameterizations")
    for p, value in zip(parameters, coefficients):
        if not np.isfinite(value) or abs(value) > p.limit + 1e-14:
            raise ValueError("Invalid coefficient")
        mask = p.mask(grid)
        coverage += mask
        change += value * mask
    if np.any(coverage > 1):
        raise ValueError("Overlapping native incident control groups")
    return values * (
        np.exp(change) if modes == {"exponential"} else 1 + change
    )
