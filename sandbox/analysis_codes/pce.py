"""Orthonormal Legendre PCE for independent uniforms, without extra dependencies."""

from dataclasses import dataclass
import itertools
import json
import numpy as np
from numpy.polynomial.legendre import legval


def indices(dimension: int, order=2, q=1.0, interaction=None) -> np.ndarray:
    """Generate total/hyperbolic degree indices without exponential enumeration."""
    if dimension < 1 or order < 0 or not 0 < q <= 1:
        raise ValueError("Invalid PCE truncation")
    result = [np.zeros(dimension, dtype=int)]
    for support in range(
        1, min(dimension, order, interaction or dimension) + 1
    ):
        for subset in itertools.combinations(range(dimension), support):
            for powers in itertools.product(
                range(1, order + 1), repeat=support
            ):
                if sum(p**q for p in powers) ** (1 / q) <= order + 1e-12:
                    a = np.zeros(dimension, dtype=int)
                    a[list(subset)] = powers
                    result.append(a)
    return np.asarray(sorted(result, key=lambda a: (sum(a), tuple(a))))


def matrix(unit: np.ndarray, powers: np.ndarray) -> np.ndarray:
    """Evaluate E[Psi_alpha Psi_beta]=delta_alpha,beta on the unit cube."""
    unit = np.asarray(unit, float)
    if (
        unit.ndim != 2
        or unit.shape[1] != powers.shape[1]
        or not np.all(np.isfinite(unit))
        or np.any((unit < 0) | (unit > 1))
    ):
        raise ValueError("Inputs must lie in the unit cube")
    basis = np.ones((len(unit), len(powers)))
    for j in range(unit.shape[1]):
        for degree in np.unique(powers[:, j]):
            mask = powers[:, j] == degree
            basis[:, mask] *= (
                np.sqrt(2 * degree + 1)
                * legval(2 * unit[:, j] - 1, [0] * degree + [1])
            )[:, None]
    return basis


@dataclass
class PCE:
    """Multiresponse PCE and training-noise coefficient covariances."""

    powers: np.ndarray
    coefficients: np.ndarray
    covariance: np.ndarray
    diagnostics: dict

    def predict(self, unit: np.ndarray) -> np.ndarray:
        return matrix(unit, self.powers) @ self.coefficients

    @property
    def mean(self):
        return self.coefficients[0]

    @property
    def variance(self):
        return np.sum(self.coefficients[1:] ** 2, axis=0)

    def save(self, path) -> None:
        np.savez_compressed(
            path,
            powers=self.powers,
            coefficients=self.coefficients,
            covariance=self.covariance,
            diagnostics=np.array(json.dumps(self.diagnostics)),
        )

    @classmethod
    def load(cls, path):
        """Read saved basis, coefficients and noise diagnostics."""
        with np.load(path) as a:
            return cls(
                a["powers"],
                a["coefficients"],
                a["covariance"],
                json.loads(str(a["diagnostics"])),
            )


def fit(
    unit: np.ndarray,
    response: np.ndarray,
    errors=None,
    *,
    order=2,
    q=1.0,
    interaction=None,
    ridge=0.0,
    sparse=False,
    keep_main_effects=False,
) -> PCE:
    """Fit response bins with weighted least squares or explicit ridge.

    Dense unregularized fits require at least twice as many samples as terms.
    Optional sparse mode iteratively prunes weak terms using TRAINING data;
    final validation is never used for model selection. Coefficient covariance
    propagates supplied tally noise, not truncation error or covariance UQ.
    """
    y = np.asarray(response, float)
    if y.ndim == 1:
        y = y[:, None]
    errors = (
        np.ones_like(y)
        if errors is None
        else np.asarray(errors, float).reshape(y.shape)
    )
    if (
        not np.all(np.isfinite(y))
        or np.any(np.isnan(errors))
        or np.any(errors < 0)
    ):
        raise ValueError("Nonfinite responses or invalid errors")
    powers = indices(unit.shape[1], order, q, interaction)
    a = matrix(unit, powers)
    if ridge < 0 or (ridge == 0 and len(unit) < 2 * len(powers)):
        raise ValueError(
            f"Need >= {2*len(powers)} samples or explicit positive regularization"
        )
    coef, cov, conditions, unavailable = [], [], [], []
    for column in range(y.shape[1]):
        s = errors[:, column]
        finite = np.isfinite(s)
        if finite.sum() < (2 * len(powers) if ridge == 0 else 2):
            coef.append(np.zeros(len(powers)))
            cov.append(np.zeros((len(powers), len(powers))))
            conditions.append(float("inf"))
            unavailable.append(column)
            continue
        positive = s[(s > 0) & finite]
        floor = max(
            np.median(positive) * 0.1 if len(positive) else 1.0, 1e-300
        )
        weights = 1 / np.maximum(s, floor)
        b = a * weights[:, None]
        conditions.append(float(np.linalg.cond(b)))
        if ridge == 0 and np.linalg.matrix_rank(b) < len(powers):
            raise ValueError("Rank-deficient PCE design")
        penalty = np.eye(len(powers)) * ridge
        penalty[0, 0] = 0
        inverse = np.linalg.pinv(b.T @ b + penalty)
        estimator = inverse @ b.T * weights[None, :]
        c = estimator @ y[:, column]
        if sparse:
            keep = np.r_[
                True, abs(c[1:]) > max(np.max(abs(c[1:])) * 0.02, 1e-300)
            ]
            if keep_main_effects:
                keep |= (powers>0).sum(axis=1)<=1
            small = b[:, keep]
            inv = np.linalg.pinv(small.T @ small + penalty[np.ix_(keep, keep)])
            estimator[:] = 0
            estimator[keep] = inv @ small.T * weights[None, :]
            c = estimator @ y[:, column]
        coef.append(c)
        noise = np.where(finite, s, 0.0)
        cov.append(
            (estimator * noise[None, :]) @ (estimator * noise[None, :]).T
        )
    return PCE(
        powers,
        np.stack(coef, axis=1),
        np.stack(cov),
        dict(
            terms=len(powers),
            samples=len(unit),
            condition_numbers=conditions,
            ridge=ridge,
            sparse=sparse,
            keep_main_effects=keep_main_effects,
            unavailable_bins=unavailable,
        ),
    )


def cross_validate(
    unit, response, errors, orders=(1, 2), folds=4, seed=42, **kwargs
):
    """Choose polynomial order using training-only folds; return all scores."""
    splits = np.array_split(
        np.random.default_rng(seed).permutation(len(unit)), folds
    )
    scores = {}
    for order in orders:
        residuals = []
        try:
            for held in splits:
                train = np.setdiff1d(np.arange(len(unit)), held)
                model = fit(
                    unit[train],
                    response[train],
                    errors[train],
                    order=order,
                    **kwargs,
                )
                residuals.extend(
                    (model.predict(unit[held]) - response[held]) ** 2
                )
            scores[order] = float(np.mean(residuals))
        except ValueError:
            scores[order] = float("inf")
    selected = min(scores, key=scores.get)
    if not np.isfinite(scores[selected]):
        raise ValueError("Insufficient training points for cross-validation")
    return selected, scores
