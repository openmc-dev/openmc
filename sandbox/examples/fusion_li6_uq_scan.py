#!/usr/bin/env python3
"""Native 14 MeV Li6 beam scan, spatial UQ maps and total-rate PCE.

Run from the repository root with OPENMC_CROSS_SECTIONS configured:
    OMP_NUM_THREADS=2 MPLBACKEND=Agg python -m sandbox.examples.fusion_li6_uq_scan

The independent uniform coefficient ranges are an explicit sensitivity prior,
NOT an evaluated nuclear-data covariance model. Each expensive sample uses
native OpenMC XS/distribution mutation, transport, tallies and exact rollback.
PCE approximates the log of total (n,t) rate; its probability distribution is
computed from the assumed input prior, not from repeated tally noise alone.
"""

import argparse
import csv
from itertools import product
import json
import hashlib
import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

# Permit both direct execution and importing this example from tests/notebooks.
if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, TwoSlopeNorm
import numpy as np
from numpy.polynomial.legendre import legval
from scipy.stats import gaussian_kde, qmc

import openmc
import openmc.lib
from sandbox.examples.fusion_li6_runtime_uq import (
    PARAMETERS,
    SPATIAL_MESH_DIMENSIONS,
    SPATIAL_MESH_LOWER,
    SPATIAL_MESH_UPPER,
    SWEEP_UPPER_BOUNDS,
    TrialData,
    make_model,
    run_trial,
    row,
    write_rows,
    save_figure,
    window_info,
)


def multi_indices(dimension, degree):
    """Total-degree polynomial indices, sorted by order then tuple."""
    return np.array(
        sorted(
            (
                a
                for a in product(range(degree + 1), repeat=dimension)
                if sum(a) <= degree
            ),
            key=lambda a: (sum(a), a),
        ),
        dtype=int,
    )


def pce_matrix(unit_inputs, indices):
    """Orthonormal Legendre products for independent uniform inputs on [0,1]."""
    # sqrt(2*n+1)*P_n(2*u-1) has unit second moment for u~Uniform(0,1).
    u = np.asarray(unit_inputs, dtype=float)
    matrix = np.ones((len(u), len(indices)))
    for j in range(u.shape[1]):
        for order in np.unique(indices[:, j]):
            polynomial = np.zeros(order + 1)
            polynomial[-1] = 1.0
            values = np.sqrt(2 * order + 1) * legval(
                2 * u[:, j] - 1, polynomial
            )
            matrix[:, indices[:, j] == order] *= values[:, None]
    return matrix


def fit_pce(unit_inputs, rates, standard_errors, degree):
    """Weighted log-rate fit: preserve positive predictions and use tally errors."""
    indices = multi_indices(unit_inputs.shape[1], degree)
    if len(unit_inputs) < 2 * len(indices):
        raise ValueError(
            f"Use at least {2*len(indices)} training samples for "
            f"{len(indices)} PCE coefficients"
        )
    if np.any(rates <= 0) or not np.all(np.isfinite(rates)):
        raise ValueError("PCE requires finite positive total production rates")
    relative_errors = standard_errors / rates
    if not np.all(np.isfinite(relative_errors)) or np.any(
        relative_errors <= 0
    ):
        raise ValueError(
            "Total-rate standard errors must be finite and positive"
        )
    matrix = pce_matrix(unit_inputs, indices)
    # The delta-method approximation gives sigma(log rate) ~= sigma(rate)/rate.
    # Bound excessive weights from accidentally very small batch error estimates.
    sigma = np.maximum(relative_errors, 0.1 * np.median(relative_errors))
    weighted = matrix / sigma[:, None]
    coefficients, _, rank, _ = np.linalg.lstsq(
        weighted, np.log(rates) / sigma, rcond=None
    )
    if rank != len(indices):
        raise RuntimeError("Rank-deficient PCE experimental design")
    # Coefficient covariance accounts for training tally noise ONLY; truncation
    # error is assessed separately on genuinely held-out transport samples.
    covariance = np.linalg.inv(weighted.T @ weighted)
    return indices, coefficients, covariance


def predict(unit_inputs, indices, coefficients):
    return np.exp(pce_matrix(unit_inputs, indices) @ coefficients)


def field_statistics(means, errors):
    """Separate sampled field variation from independent transport estimator noise."""
    average = np.mean(means, axis=0)
    observed_variance = np.var(means, axis=0, ddof=1)
    # Some native zero-score bins report infinite errors. Do not fabricate a
    # noise correction there: mark those voxels unavailable in corrected maps.
    available = np.all(np.isfinite(errors), axis=0)
    noise_variance = np.full(average.shape, np.nan)
    noise_variance[available] = np.mean(errors[:, available] ** 2, axis=0)
    corrected = np.maximum(observed_variance - noise_variance, 0.0)
    return (
        average,
        np.sqrt(corrected),
        np.sqrt(noise_variance),
        observed_variance,
    )


def spatial_images(
    directory,
    means,
    errors,
    baseline,
    totals,
    depth_means,
    depth_errors,
    baseline_depth,
):
    """Plot sampled 3D spatial variability; no surrogate field is fabricated."""
    mean, std, noise, observed_variance = field_statistics(means, errors)
    middle = SPATIAL_MESH_DIMENSIONS[1] // 2
    edges = [
        np.linspace(lo, hi, n + 1)
        for lo, hi, n in zip(
            SPATIAL_MESH_LOWER, SPATIAL_MESH_UPPER, SPATIAL_MESH_DIMENSIONS
        )
    ]
    # Save the entire field statistics; sample meshes are saved separately.
    np.savez_compressed(
        directory / "spatial_statistics.npz",
        mean=mean,
        parameter_std_estimate=std,
        transport_error_rms=noise,
        observed_variance=observed_variance,
        x_edges_cm=edges[0],
        y_edges_cm=edges[1],
        z_edges_cm=edges[2],
    )
    baseline_slice = baseline[:, middle, :]
    mean_slice = mean[:, middle, :]
    delta = mean_slice - baseline_slice
    maximum = max(float(np.max(abs(delta))), np.finfo(float).tiny)
    positive = np.r_[
        baseline_slice[baseline_slice > 0], mean_slice[mean_slice > 0]
    ]
    rate_norm = LogNorm(vmin=positive.min(), vmax=positive.max())

    def panels(filename, arrays, titles, norms, cmaps, labels):
        fig, axes = plt.subplots(
            1, len(arrays), figsize=(5 * len(arrays), 5), layout="constrained"
        )
        for ax, values, title, norm, cmap, label in zip(
            np.atleast_1d(axes), arrays, titles, norms, cmaps, labels
        ):
            if isinstance(norm, LogNorm):
                values = np.ma.masked_less_equal(values, 0.0)
            picture = ax.pcolormesh(
                edges[0],
                edges[2],
                np.ma.masked_invalid(values).T,
                cmap=cmap,
                norm=norm,
                shading="flat",
            )
            ax.add_patch(
                plt.Circle((0, 0), 20, fill=False, color="grey", linewidth=0.8)
            )
            ax.annotate(
                "",
                xy=(0, -16),
                xytext=(0, -20),
                arrowprops=dict(arrowstyle="->", color="black"),
            )
            ax.set(
                title=title, xlabel="x [cm]", ylabel="z [cm]; beam along +z"
            )
            ax.set_aspect("equal")
            fig.colorbar(picture, ax=ax, label=label, shrink=0.8)
        fig.suptitle(
            "14 MeV Li6 beam: central y=0 finite-thickness voxel slice"
        )
        fig.savefig(directory / filename, dpi=180)
        plt.close(fig)

    panels(
        "beam_mean.png",
        [baseline_slice, mean_slice, delta],
        [
            "Unperturbed baseline",
            "Mean across input prior",
            "Mean minus baseline",
        ],
        [
            rate_norm,
            rate_norm,
            TwoSlopeNorm(vmin=-maximum, vcenter=0, vmax=maximum),
        ],
        ["viridis", "viridis", "RdBu_r"],
        ["(n,t) / voxel / source"] * 3,
    )
    std_slice, noise_slice = std[:, middle, :], noise[:, middle, :]
    cv = np.divide(
        std_slice,
        mean_slice,
        out=np.full_like(mean_slice, np.nan),
        where=mean_slice > 0,
    )
    # Mask very weak-production cells, where relative uncertainty is dominated
    # by finite histories. Grey/uncoloured cells also denote unavailable errors.
    cv[mean_slice < mean_slice.max() * 1.0e-4] = np.nan
    panels(
        "beam_variation.png",
        [std_slice, 100 * cv, noise_slice],
        [
            "Parameter SD (transport-noise corrected)",
            "Parameter coefficient of variation",
            "RMS transport standard error",
        ],
        [None, None, None],
        ["magma"] * 3,
        ["(n,t) / voxel / source", "[%]", "(n,t) / voxel / source"],
    )
    selected = [
        int(np.argmin(totals)),
        int(np.argsort(totals)[len(totals) // 2]),
        int(np.argmax(totals)),
    ]
    selected_slices = [means[i, :, middle, :] for i in selected]
    positive = np.concatenate([a[a > 0] for a in selected_slices])
    shared = LogNorm(vmin=positive.min(), vmax=positive.max())
    panels(
        "beam_sample_extremes.png",
        selected_slices,
        [f"Sample {i}: total={totals[i]:.5f}" for i in selected],
        [shared] * 3,
        ["viridis"] * 3,
        ["(n,t) / voxel / source"] * 3,
    )
    average, parameter_std, _, _ = field_statistics(depth_means, depth_errors)
    depth = np.linspace(1, 39, 20)
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(depth, baseline_depth, label="Unperturbed baseline")
    ax.plot(depth, average, label="Mean across input prior")
    ax.fill_between(
        depth,
        average - parameter_std,
        average + parameter_std,
        alpha=0.25,
        label="Estimated parameter SD (noise corrected)",
    )
    ax.set(
        xlabel="Depth from beam entrance [cm]", ylabel="(n,t) / slab / source"
    )
    ax.legend()
    ax.grid(alpha=0.2)
    save_figure(fig, directory / "beam_depth_variation.png")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--profile",
        choices=["symmetric10", "symmetric20", "log-tilt"],
        default="symmetric10",
    )
    parser.add_argument(
        "--baseline-particles",
        type=int,
        default=400000,
        help="particles per batch for the high-statistics common baseline",
    )
    parser.add_argument("--train-samples", type=int, default=64)
    parser.add_argument("--validation-samples", type=int, default=32)
    parser.add_argument("--particles", type=int, default=40000)
    parser.add_argument("--batches", type=int, default=20)
    parser.add_argument("--degree", type=int, default=None)
    parser.add_argument("--surrogate-draws", type=int, default=100000)
    parser.add_argument("--seed", type=int, default=20261014)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("sandbox/examples/fusion_li6_uq_symmetric10_output"),
    )
    parser.add_argument(
        "--lower-bounds",
        type=float,
        nargs=5,
        default=None,
        metavar=("A", "B1", "B2", "C1", "C2"),
    )
    parser.add_argument(
        "--upper-bounds",
        type=float,
        nargs=5,
        default=None,
        metavar=("A", "B1", "B2", "C1", "C2"),
    )
    parser.add_argument(
        "--refit-only",
        type=Path,
        help="Refit saved native scan data without rerunning transport",
    )
    args = parser.parse_args()
    cached = None
    if args.refit_only is not None:
        if args.lower_bounds is not None or args.upper_bounds is not None:
            parser.error(
                "Cannot change input prior bounds when refitting cached transport"
            )
        args.output = args.refit_only
        cached = json.loads((args.output / "metadata.json").read_text())
        args.profile = cached.get("perturbation_profile", "log-tilt")
        args.train_samples = cached["scan"]["training"]
        args.validation_samples = cached["scan"]["validation"]
        args.particles = cached["scan"]["particles_per_batch"]
        args.batches = cached["scan"]["batches"]
        args.seed = cached["scan"]["design_seed"]
        args.baseline_particles = (
            cached["scan"].get(
                "baseline_histories", args.particles * args.batches
            )
            // args.batches
        )
        args.lower_bounds = cached["prior"]["lower_bounds"]
        args.upper_bounds = cached["prior"]["upper_bounds"]
    if args.degree is None:
        args.degree = (
            1 if args.profile in ("symmetric10", "symmetric20") else 2
        )
    terms = len(multi_indices(5, args.degree)) if args.degree >= 0 else 0
    if (
        args.degree < 1
        or args.train_samples < 2 * terms
        or args.validation_samples < 8
        or args.particles < 1
        or args.batches < 2
        or args.surrogate_draws < 1000
    ):
        parser.error(
            "Require degree>=1, training samples>=2*PCE terms, validation>=8, "
            "positive particles, batches>=2 and surrogate draws>=1000"
        )
    args.output.mkdir(parents=True, exist_ok=True)
    from sandbox.examples.li6_symmetric_perturbations import (
        SymmetricTrialData,
        paired_design,
        paired_trial,
        PARAMETERS as FRACTION_NAMES,
    )

    symmetric = args.profile in ("symmetric10", "symmetric20")
    fraction_limit = 0.2 if args.profile == "symmetric20" else 0.1
    parameter_names = FRACTION_NAMES if symmetric else PARAMETERS
    lower = np.array(
        args.lower_bounds
        if args.lower_bounds is not None
        else ([-fraction_limit] * 5 if symmetric else [0.0] * 5)
    )
    bounds = np.array(
        args.upper_bounds
        if args.upper_bounds is not None
        else ([fraction_limit] * 5 if symmetric else SWEEP_UPPER_BOUNDS)
    )
    if symmetric and (
        not np.array_equal(lower, -bounds)
        or np.any(abs(bounds) > fraction_limit)
        or args.train_samples % 2
        or args.validation_samples % 2
    ):
        parser.error(
            f"{args.profile} requires symmetric bounds within +/-{fraction_limit} and even sample counts"
        )
    if (
        not np.all(np.isfinite(lower))
        or not np.all(np.isfinite(bounds))
        or np.any(bounds <= lower)
    ):
        parser.error(
            "Every upper bound must be finite and greater than its finite lower bound"
        )
    # Independent designs and unique transport seeds avoid training/validation
    # leakage and allow a real predictive check of the surrogate.
    train = (
        paired_design(args.train_samples, args.seed)
        if symmetric
        else qmc.LatinHypercube(5, seed=args.seed).random(args.train_samples)
    )
    validation = (
        paired_design(args.validation_samples, args.seed + 1)
        if symmetric
        else qmc.LatinHypercube(5, seed=args.seed + 1).random(
            args.validation_samples
        )
    )
    units = np.vstack((train, validation))
    inputs = lower + units * (bounds - lower)
    if cached is None:
        totals = []
        total_errors = []
        meshes = []
        mesh_errors = []
        depths = []
        depth_errors = []
        records = []
        with TemporaryDirectory(prefix="li6-full-uq-") as temporary:
            model, ids = make_model(args.particles, args.batches, temporary)
            try:
                model.init_lib(directory=temporary, output=False, threads=1)
                data = (
                    SymmetricTrialData(fraction_limit)
                    if symmetric
                    else TrialData()
                )
                if symmetric:
                    openmc.lib.settings.particles = args.baseline_particles
                baseline_result = run_trial(data, ids, [0.0] * 5, 80000)
                openmc.lib.settings.particles = args.particles
                if symmetric:
                    print(
                        f"Symmetric fractional prior: XS multiplier {1-fraction_limit:g}..{1+fraction_limit:g}; two shape modes share a +/-{100*fraction_limit:g}% pointwise PDF budget. Matched-history differences; baseline SE is separate.",
                        flush=True,
                    )
                baseline = baseline_result["tritium_mesh"][0].reshape(
                    SPATIAL_MESH_DIMENSIONS, order="F"
                )
                for i, parameters in enumerate(inputs):
                    seed = 100000 + i
                    label = "train" if i < args.train_samples else "validation"
                    try:
                        if symmetric:
                            result, paired_details = paired_trial(
                                data, ids, parameters, seed, baseline_result
                            )
                        else:
                            result = run_trial(data, ids, parameters, seed)
                            paired_details = {}
                    except Exception as error:
                        failure = row(
                            f"{label}_{i}",
                            parameters,
                            seed,
                            args.particles,
                            args.batches,
                            error=str(error),
                        )
                        failure["split"] = label
                        # Save evidence and stop; a failed sample must not silently
                        # disappear or be replaced by a fabricated surrogate target.
                        (args.output / "failure.json").write_text(
                            json.dumps(failure, indent=2)
                        )
                        raise
                    totals.append(float(result["tritium_rate"][0][0]))
                    total_errors.append(float(result["tritium_rate"][1][0]))
                    meshes.append(
                        result["tritium_mesh"][0].reshape(
                            SPATIAL_MESH_DIMENSIONS, order="F"
                        )
                    )
                    mesh_errors.append(
                        result["tritium_mesh"][1].reshape(
                            SPATIAL_MESH_DIMENSIONS, order="F"
                        )
                    )
                    depths.append(result["tritium_depth"][0])
                    depth_errors.append(result["tritium_depth"][1])
                    record = row(
                        f"{label}_{i}",
                        parameters,
                        seed,
                        args.particles,
                        args.batches,
                        result,
                    )
                    record.update(
                        split=label,
                        energy_mt=data.product["mt"],
                        energy_product=data.product["product_index"],
                        energy_distribution=data.descriptor[
                            "distribution_index"
                        ],
                    )
                    if symmetric:
                        for old, new in zip(PARAMETERS, parameter_names):
                            record[new] = record.pop(old)
                        record.update(paired_details)
                    records.append(record)
                    if (i + 1) % 8 == 0 or i + 1 == len(inputs):
                        print(
                            f"{i+1}/{len(inputs)} native runs: "
                            f"total (n,t)={totals[-1]:.6g} +/- {total_errors[-1]:.2g}",
                            flush=True,
                        )
                        write_rows(args.output / "scan.csv", records)
                if symmetric:
                    openmc.lib.settings.particles = args.baseline_particles
                restored = run_trial(data, ids, [0.0] * 5, 80000)
                for name in ids:
                    for a, b in zip(baseline_result[name], restored[name]):
                        np.testing.assert_array_equal(a, b)
                selection = dict(
                    xs_mt=2,
                    angle_mt=2,
                    energy_mt=data.product["mt"],
                    energy_product=data.product["product_index"],
                    energy_distribution=data.descriptor["distribution_index"],
                )
                windows = dict(
                    xs=window_info(data.grid),
                    angle=window_info(data.angle["incident_energy"]),
                    energy=window_info(data.energy["incident_energy"]),
                )
            finally:
                model.finalize_lib()
        totals = np.array(totals)
        total_errors = np.array(total_errors)
        means = np.stack(meshes)
        errors = np.stack(mesh_errors)
        np.savez_compressed(
            args.output / "scan_meshes.npz",
            parameters=inputs,
            unit_inputs=units,
            mean=means,
            std_error=errors,
            total=totals,
            total_std_error=total_errors,
            baseline_mean=baseline,
            dimensions=SPATIAL_MESH_DIMENSIONS,
            lower_left_cm=SPATIAL_MESH_LOWER,
            upper_right_cm=SPATIAL_MESH_UPPER,
            split=np.array([r["split"] for r in records]),
        )
        spatial_images(
            args.output,
            means,
            errors,
            baseline,
            totals,
            np.stack(depths),
            np.stack(depth_errors),
            baseline_result["tritium_depth"][0],
        )
    else:
        with np.load(args.output / "scan_meshes.npz") as saved:
            units = saved["unit_inputs"].copy()
            inputs = saved["parameters"].copy()
            totals = saved["total"].copy()
            total_errors = saved["total_std_error"].copy()
        train = units[: args.train_samples]
        validation = units[args.train_samples :]
        baseline_result = {
            "tritium_rate": (
                np.array(
                    [cached["pce"]["predicted_distribution"]["baseline_mean"]]
                ),
                np.array([cached.get("common_baseline_standard_error", 0.0)]),
            )
        }
        selection = cached["selection"]
        windows = cached["incident_windows"]
        print(
            "Refitting cached native scan; no OpenMC initialization or transport.",
            flush=True,
        )
    n = args.train_samples
    indices, coefficients, covariance = fit_pce(
        train, totals[:n], total_errors[:n], args.degree
    )
    predicted = predict(validation, indices, coefficients)
    residual = predicted - totals[n:]
    r2 = 1 - np.sum(residual**2) / np.sum(
        (totals[n:] - np.mean(totals[n:])) ** 2
    )
    rmse = float(np.sqrt(np.mean(residual**2)))
    validation_matrix = pce_matrix(validation, indices)
    prediction_se = predicted * np.sqrt(
        np.maximum(
            0.0,
            np.einsum(
                "ij,jk,ik->i", validation_matrix, covariance, validation_matrix
            ),
        )
    )
    normalized_residual = residual / np.sqrt(
        total_errors[n:] ** 2 + prediction_se**2
    )
    parameter_spread = np.std(
        predict(
            np.random.default_rng(args.seed + 3).random((20000, 5)),
            indices,
            coefficients,
        )
    )
    fit_response_noise = float(
        np.mean(totals[:n]) * np.sqrt(np.trace(covariance[1:, 1:]))
    )
    noise_aware_accept = bool(
        np.mean(normalized_residual**2) <= 2.0
        and parameter_spread > 3 * fit_response_noise
        and rmse / np.mean(totals[n:]) <= 0.005
    )
    validation_metrics = dict(
        r2=float(r2),
        rmse=rmse,
        relative_rmse=rmse / float(np.mean(totals[n:])),
        rms_transport_standard_error=float(
            np.sqrt(np.mean(total_errors[n:] ** 2))
        ),
        samples=args.validation_samples,
        accepted=(
            noise_aware_accept
            if symmetric
            else bool(r2 >= 0.9 and rmse / np.mean(totals[n:]) <= 0.02)
        ),
        method="paired-noise-aware" if symmetric else "raw-r2",
        mean_squared_standardized_residual=float(
            np.mean(normalized_residual**2)
        ),
        parameter_spread_to_fit_noise=float(
            parameter_spread / fit_response_noise
        ),
    )
    np.savez(
        args.output / "pce_model.npz",
        multi_indices=indices,
        coefficients=coefficients,
        coefficient_covariance_from_tally_noise=covariance,
        lower_bounds=lower,
        upper_bounds=bounds,
        log_response=True,
        common_baseline_standard_error=(
            baseline_result["tritium_rate"][1][0] if symmetric else 0.0
        ),
    )
    write_rows(
        args.output / "pce_validation.csv",
        [
            dict(
                transport_mean=float(y),
                transport_std_error=float(e),
                pce_prediction=float(p),
            )
            for y, e, p in zip(totals[n:], total_errors[n:], predicted)
        ],
    )
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.errorbar(
        totals[n:],
        predicted,
        xerr=total_errors[n:],
        fmt="o",
        alpha=0.75,
        label="Held-out transport (horizontal bars: 1 SE)",
    )
    limits = [
        min(totals[n:].min(), predicted.min()),
        max(totals[n:].max(), predicted.max()),
    ]
    ax.plot(limits, limits, "k--")
    ax.legend(fontsize=8)
    ax.set(
        xlabel="Native transport (n,t) / source",
        ylabel="PCE prediction",
        title=f"Independent validation: R²={r2:.3f}, RMSE={rmse:.3g}",
    )
    save_figure(fig, args.output / "pce_validation.png")
    # Propagate the assumed parameter prior through the fitted response, not
    # through an independent particle sampler or nuclear-data processing engine.
    draw_units = np.random.default_rng(args.seed + 2).random(
        (args.surrogate_draws, 5)
    )
    draws = predict(draw_units, indices, coefficients)
    pdf_summary = dict(
        mean=float(draws.mean()),
        standard_deviation=float(draws.std(ddof=1)),
        quantiles_5_50_95=np.quantile(draws, [0.05, 0.5, 0.95]).tolist(),
        baseline_mean=float(baseline_result["tritium_rate"][0][0]),
        surrogate_draws=args.surrogate_draws,
    )
    # Never publish a density as validated UQ when held-out transport disagrees.
    # Save diagnostics/fields first so an unsuccessful fit remains reviewable.
    library = Path(os.environ["OPENMC_CROSS_SECTIONS"])
    metadata = dict(
        openmc_version=openmc.__version__,
        git_commit=subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True
        ).strip(),
        git_dirty=bool(
            subprocess.check_output(
                ["git", "status", "--porcelain"], text=True
            ).strip()
        ),
        library=str(library),
        library_index_sha256=hashlib.sha256(library.read_bytes()).hexdigest(),
        parameters=list(parameter_names),
        perturbation_profile=args.profile,
        refit_from_cached_native_scan=bool(cached is not None),
        common_baseline_standard_error=(
            float(baseline_result["tritium_rate"][1][0]) if symmetric else 0.0
        ),
        prior=dict(
            type=(
                "independent uniform fractional prior"
                if symmetric
                else "independent uniform sensitivity prior"
            ),
            design=(
                "antithetic Latin hypercube with exact +/- input pairs"
                if symmetric
                else "Latin hypercube"
            ),
            pointwise_pdf_budget=fraction_limit if symmetric else None,
            lower_bounds=lower.tolist(),
            upper_bounds=bounds.tolist(),
            evaluated_covariance=False,
        ),
        selection=selection,
        incident_windows=windows,
        basis=dict(
            xs=(
                f"sigma*(1+delta), delta uniform(-{fraction_limit},{fraction_limit})"
                if symmetric
                else "log scale at native 10–14 MeV grid points"
            ),
            energy=(
                "Native-mass-centered, bounded z and z² modes; 1+(b1*g1+b2*g2)/2"
                if symmetric
                else "z=E_out/max(table E_out); b1*z+b2*z²"
            ),
            angle=(
                "Native-mass-centered, bounded P1 and P2 modes; 1+(c1*g1+c2*g2)/2"
                if symmetric
                else "c1*mu+c2*(3mu²-1)/2"
            ),
        ),
        source=dict(
            energy_eV=14.0e6,
            position_cm=[0, 0, -20 + 1.0e-6],
            direction=[0, 0, 1],
        ),
        geometry=dict(radius_cm=20, density_g_cm3=0.46),
        scan=dict(
            training=args.train_samples,
            validation=args.validation_samples,
            histories_per_sample=args.particles * args.batches,
            particles_per_batch=args.particles,
            batches=args.batches,
            design_seed=args.seed,
            transport_seeds_start=100000,
            baseline_seed=80000,
            transport_threads=1,
            exact_restoration=True,
            baseline_histories=(
                args.baseline_particles * args.batches
                if symmetric
                else args.particles * args.batches
            ),
            matched_baseline_runs=symmetric,
            native_scan_transport_runs=(
                2 * len(inputs) if symmetric else len(inputs)
            ),
            native_scan_histories=(2 if symmetric else 1)
            * len(inputs)
            * args.particles
            * args.batches,
            estimator=(
                "common high-stat baseline + paired perturbed-minus-control mean"
                if symmetric
                else "native tally mean"
            ),
            standard_errors=(
                "paired batch-difference SE, excludes separately reported common baseline SE"
                if symmetric
                else "native tally SE"
            ),
        ),
        pce=dict(
            degree=args.degree,
            terms=len(indices),
            response="log total MT105 rate/source",
            basis="orthonormal Legendre products",
            weighted_by="relative tally standard errors",
            validation=validation_metrics,
            predicted_distribution=pdf_summary,
        ),
        spatial=dict(
            units="(n,t) reactions / voxel / source",
            dimensions=SPATIAL_MESH_DIMENSIONS,
            order="[sample,x,y,z], native x-fastest/Fortran reshape",
            variation="sample variance minus average transport variance, clipped at zero",
            unavailable="voxels with nonfinite native errors retain NaN in corrected statistics",
            caveat="sampled spatial extremes and quantiles include transport estimator noise",
        ),
        limitations=[
            "Prior is user-assumed, not evaluated nuclear-data uncertainties.",
            "PCE is a validated approximation for total rate, not a field or native sampler replacement.",
            "Spatial variation has finite-design error; corrected variance can be unresolved in noisy voxels.",
            "PDF excludes transport estimator noise; PCE truncation and fit error remain as validated.",
            "Symmetric profile PDF conditions on a common baseline estimate; its SE is reported separately. Raw perturbed/control totals are retained in scan.csv.",
        ],
    )
    (args.output / "metadata.json").write_text(json.dumps(metadata, indent=2))
    if not validation_metrics["accepted"]:
        raise RuntimeError(
            f"PCE failed held-out validation: {validation_metrics}; "
            "increase histories/design size or polynomial degree before using its PDF"
        )
    padding = 0.1 * (draws.max() - draws.min())
    grid = np.linspace(draws.min() - padding, draws.max() + padding, 500)
    density = gaussian_kde(draws)(grid)
    cdf = np.searchsorted(np.sort(draws), grid, side="right") / len(draws)
    write_rows(
        args.output / "total_nt_distribution.csv",
        [
            dict(rate_per_source=float(x), pdf=float(p), cdf=float(c))
            for x, p, c in zip(grid, density, cdf)
        ],
    )
    np.savez_compressed(
        args.output / "pce_distribution_samples.npz", rate_per_source=draws
    )
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    axes[0].plot(grid, density, label="Validated log-rate PCE + assumed prior")
    axes[0].hist(
        totals,
        bins=16,
        density=True,
        alpha=0.25,
        label="Transport response estimates (include estimator noise)",
    )
    axes[0].axvline(
        pdf_summary["baseline_mean"],
        color="k",
        linestyle="--",
        label="Unperturbed baseline",
    )
    axes[0].set(
        xlabel="Total Li6 (n,t) reactions / source neutron",
        ylabel="Probability density",
    )
    axes[0].legend(fontsize=8)
    axes[1].plot(grid, cdf)
    axes[1].set(
        xlabel="Total Li6 (n,t) reactions / source neutron",
        ylabel="Cumulative probability",
        ylim=(0, 1),
    )
    fig.suptitle(
        "Parameter-induced production distribution (assumed independent uniform prior)"
    )
    save_figure(fig, args.output / "total_nt_distribution.png")
    print("PCE validation:", validation_metrics, flush=True)
    print("Predicted (n,t) distribution:", pdf_summary, flush=True)
    print("Output:", args.output)


if __name__ == "__main__":
    main()
