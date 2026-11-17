"""Symmetric fractional changes for the Li6 runtime UQ example.

XS multiplier is 1+delta, not exp(delta). Shape changes are linear in zero-mean
coefficients and centered with the native interpolation rule before mutation.
Two bounded shape modes share the configured fractional budget, so their combination remains
within that budget pointwise (apart from native normalization roundoff).
"""

from copy import deepcopy

import numpy as np
import openmc.lib
from sandbox.examples.fusion_li6_runtime_uq import (
    TrialData,
    window_mask,
    XS_MT,
    energy_basis,
    angle_basis,
)
from openmc.lib.perturbation import _areas

PARAMETERS = (
    "xs_fraction",
    "energy_fraction_1",
    "energy_fraction_2",
    "angle_fraction_1",
    "angle_fraction_2",
)


def symmetric_shape(law, fractions, basis, limit=0.1):
    fractions = np.asarray(fractions, dtype=float)
    if (
        fractions.shape != (2,)
        or not np.all(np.isfinite(fractions))
        or np.any(abs(fractions) > limit + 1.0e-14)
    ):
        raise ValueError(
            f"Two finite fractional shape parameters within +/-{limit} required"
        )
    result = deepcopy(law)
    if np.all(fractions == 0):
        return result
    key = "energy" if law["type"] == "continuous-tabular" else "x"
    for i in np.flatnonzero(window_mask(law["incident_energy"])):
        table = result["tables"][i]
        n = table.get("n_discrete", 0)
        x = table[key][n:]
        p = table["pdf"][n:].copy()
        if len(x) < 2:
            continue
        interp = table["interpolation"]
        mass = np.sum(_areas(x, p, interp))
        target = 1.0 - sum(table["pdf"][:n])
        functions = np.asarray(basis(table[key], law["incident_energy"][i]))[
            :, n:
        ]
        modes = []
        for function in functions:
            # Integrate the native interpolant of p*h, not a different smooth
            # product between its tabulated points. This makes its mass zero.
            center = np.sum(_areas(x, p * function, interp)) / mass
            centered = function - center
            maximum = np.max(abs(centered))
            modes.append(
                centered / maximum if maximum > 0 else np.zeros_like(centered)
            )
        change = 0.5 * fractions @ np.array(modes)
        table["pdf"][n:] = p * (1.0 + change) * (target / mass)
        table["cdf"][:n] = np.cumsum(table["pdf"][:n])
        table["cdf"][n:] = (
            sum(table["pdf"][:n])
            + np.r_[0.0, np.cumsum(_areas(x, table["pdf"][n:], interp))]
        )
    return result


class SymmetricTrialData(TrialData):
    def __init__(self, fraction_limit=0.1):
        if not np.isfinite(fraction_limit) or not 0 < fraction_limit < 1:
            raise ValueError(
                "Fractional budget must be finite and between zero and one"
            )
        self.fraction_limit = fraction_limit
        super().__init__()

    def candidates(self, parameters):
        parameters = np.asarray(parameters, dtype=float)
        if (
            parameters.shape != (5,)
            or not np.all(np.isfinite(parameters))
            or np.any(abs(parameters) > self.fraction_limit + 1.0e-14)
        ):
            raise ValueError(
                f"Five fractional parameters within +/-{self.fraction_limit} required"
            )
        a, b1, b2, c1, c2 = parameters
        result = dict(xs=[], energy=[], angle=[])
        if a:
            result["xs"] = [
                dict(
                    nuclide="Li6",
                    mt=XS_MT,
                    values=self.full_xs * (1.0 + a * window_mask(self.grid)),
                )
            ]
        if b1 or b2:
            result["energy"] = [
                dict(
                    nuclide="Li6",
                    mt=self.product["mt"],
                    product=self.product["product_index"],
                    distribution=self.descriptor["distribution_index"],
                    law=symmetric_shape(
                        self.energy,
                        [b1, b2],
                        energy_basis,
                        self.fraction_limit,
                    ),
                )
            ]
        if c1 or c2:
            result["angle"] = [
                dict(
                    nuclide="Li6",
                    mt=2,
                    law=symmetric_shape(
                        self.angle, [c1, c2], angle_basis, self.fraction_limit
                    ),
                )
            ]
        return result


def paired_design(count, seed):
    """Antithetic Latin hypercube: each unit vector u has a partner 1-u."""
    from scipy.stats import qmc

    if count % 2:
        raise ValueError("An even sample count is required for balanced pairs")
    half = qmc.LatinHypercube(5, seed=seed).random(count // 2)
    return np.vstack((half, 1.0 - half))


def batch_trial(data, ids, parameters, seed):
    """Collect batch estimates so paired transport error includes covariance."""
    data.assert_restored()
    openmc.lib.hard_reset()
    openmc.lib.settings.seed = seed
    estimates = {name: [] for name in ids}
    previous = {name: 0.0 for name in ids}
    with openmc.lib.nuclear_data_perturbation(**data.candidates(parameters)):
        try:
            for _ in range(openmc.lib.settings.get_batches()):
                openmc.lib.next_batch()
                for name, tid in ids.items():
                    tally = openmc.lib.tallies[tid]
                    accumulated = tally.results[:, :, 1].ravel().copy()
                    estimates[name].append(accumulated - previous[name])
                    previous[name] = accumulated
        finally:
            openmc.lib.simulation_finalize()
    data.assert_restored()
    return {name: np.stack(values) for name, values in estimates.items()}


def paired_trial(data, ids, parameters, seed, anchor):
    control = batch_trial(data, ids, [0.0] * 5, seed)
    changed = batch_trial(data, ids, parameters, seed)
    result = {}
    for name in ids:
        difference = changed[name] - control[name]
        result[name] = (
            anchor[name][0] + difference.mean(axis=0),
            difference.std(axis=0, ddof=1) / np.sqrt(len(difference)),
        )
    details = dict(
        native_perturbed_total=float(changed["tritium_rate"].mean()),
        native_control_total=float(control["tritium_rate"].mean()),
        paired_delta=float(
            (changed["tritium_rate"] - control["tritium_rate"]).mean()
        ),
        paired_delta_standard_error=float(result["tritium_rate"][1][0]),
        common_baseline_mean=float(anchor["tritium_rate"][0][0]),
        common_baseline_standard_error=float(anchor["tritium_rate"][1][0]),
    )
    return result, details
