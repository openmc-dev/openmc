#!/usr/bin/env python3
"""Increase Li6 UQ sampling precision using a validated saved PCE.

No OpenMC initialization, transport, or nuclear-data mutation is performed.
The input prior and PCE remain fixed; more draws improve propagation statistics,
not the surrogate's accuracy or its native training-data statistics.
"""

import argparse
import json
from pathlib import Path
import sys

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde
from sandbox.examples.fusion_li6_uq_scan import pce_matrix, save_figure


def propagate(directory, draws, seed, output):
    metadata = json.loads((directory / "metadata.json").read_text())
    if not metadata["pce"]["validation"]["accepted"]:
        raise ValueError(
            "PCE has not passed independent native transport validation"
        )
    with np.load(directory / "pce_model.npz") as model:
        indices = model["multi_indices"].copy()
        coefficients = model["coefficients"].copy()
        covariance = model["coefficient_covariance_from_tally_noise"].copy()
    # Use bounded blocks to avoid allocating a million-by-21 matrix at once.
    rng = np.random.default_rng(seed)
    values = np.empty(draws)
    gradient = np.zeros(len(coefficients))
    for start in range(0, draws, 20000):
        end = min(start + 20000, draws)
        basis = pce_matrix(rng.random((end - start, 5)), indices)
        response = np.exp(basis @ coefficients)
        values[start:end] = response
        # dE[exp(PCE)]/dcoef = E[exp(PCE)*basis]. Together with fitted
        # coefficient covariance this estimates native training-noise effects
        # on the predicted mean. It excludes PCE truncation/model bias.
        gradient += response @ basis
    gradient /= draws
    anchor_error = float(metadata.get("common_baseline_standard_error", 0.0))
    mean = float(values.mean())
    sd = float(values.std(ddof=1))
    quantiles = np.quantile(values, [0.05, 0.5, 0.95])
    kernel = gaussian_kde(values)
    quantile_errors = np.sqrt(
        np.array([0.05, 0.5, 0.95]) * (1 - np.array([0.05, 0.5, 0.95])) / draws
    ) / kernel(quantiles)
    summary = dict(
        draws=draws,
        seed=seed,
        mean=mean,
        parameter_standard_deviation=sd,
        mean_sampling_standard_error=sd / np.sqrt(draws),
        quantiles_5_50_95=quantiles.tolist(),
        quantile_sampling_standard_errors=quantile_errors.tolist(),
        mean_fit_standard_error_from_training_tally_noise=float(
            np.sqrt(max(0.0, gradient @ covariance @ gradient))
        ),
        common_baseline_standard_error=anchor_error,
        mean_fit_standard_error_including_common_baseline=float(
            np.sqrt(
                max(0.0, gradient @ covariance @ gradient) + anchor_error**2
            )
        ),
        native_validation=metadata["pce"]["validation"],
        assumptions="Same independent uniform prior and validated saved log-rate PCE.",
        limitation="Sampling SEs condition on the fitted PCE. More draws do not remove fit noise or truncation bias.",
    )
    output.mkdir(parents=True, exist_ok=True)
    (output / "summary.json").write_text(json.dumps(summary, indent=2))
    np.savez_compressed(output / "samples.npz", rate_per_source=values)
    padding = 0.1 * (values.max() - values.min())
    grid = np.linspace(values.min() - padding, values.max() + padding, 500)
    density = kernel(grid)
    cumulative = np.searchsorted(np.sort(values), grid, side="right") / draws
    np.savetxt(
        output / "distribution.csv",
        np.column_stack((grid, density, cumulative)),
        delimiter=",",
        header="rate_per_source,pdf,cdf",
        comments="",
    )
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    axes[0].plot(grid, density)
    axes[0].axvline(
        metadata["pce"]["predicted_distribution"]["baseline_mean"],
        color="k",
        linestyle="--",
        label="Unperturbed baseline",
    )
    axes[0].set(
        xlabel="Total Li6 (n,t) reactions/source", ylabel="Probability density"
    )
    axes[0].legend()
    axes[1].plot(grid, cumulative)
    axes[1].set(
        xlabel="Total Li6 (n,t) reactions/source",
        ylabel="Cumulative probability",
        ylim=(0, 1),
    )
    fig.suptitle(
        f"{draws:,} draws from validated PCE; no additional transport"
    )
    save_figure(fig, output / "distribution.png")
    print(json.dumps(summary, indent=2))
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--model-directory",
        type=Path,
        default=Path("sandbox/examples/fusion_li6_uq_symmetric10_output"),
    )
    parser.add_argument("--draws", type=int, default=1000000)
    parser.add_argument("--seed", type=int, default=20261016)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.draws < 1000:
        parser.error("Use at least 1000 draws")
    propagate(
        args.model_directory,
        args.draws,
        args.seed,
        args.output or args.model_directory / "pce_propagation",
    )


if __name__ == "__main__":
    main()
