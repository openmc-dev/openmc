#!/usr/bin/env python3
"""One-at-a-time native Li6 perturbations and their (n,t) transport responses.

Run from the repository root. Two shape coefficients move together; this is
one documented direction through each two-dimensional shape parameter space,
not a marginal distribution over the UQ prior. Nothing changes MT105 directly.
"""

import argparse
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import openmc.lib
from sandbox.examples.fusion_li6_runtime_uq import (
    make_model,
    run_trial,
    SPATIAL_MESH_DIMENSIONS,
)
from sandbox.examples.li6_symmetric_perturbations import (
    SymmetricTrialData,
    paired_trial,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fraction", type=float, default=0.2)
    parser.add_argument("--particles", type=int, default=40000)
    parser.add_argument("--batches", type=int, default=20)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "sandbox/examples/fusion_li6_uq_symmetric20_output/separate_effects"
        ),
    )
    args = parser.parse_args()
    if not 0 < args.fraction < 1 or args.particles < 1 or args.batches < 2:
        parser.error("Require 0<fraction<1, particles>0 and batches>=2")
    args.output.mkdir(parents=True, exist_ok=True)
    groups = {
        "cross_section": [0],
        "emission_energy": [1, 2],
        "emission_angle": [3, 4],
    }
    labels = {
        "cross_section": "Li6 MT2 elastic cross section",
        "emission_energy": "Discovered Li6 neutron energy shape",
        "emission_angle": "Li6 MT2 elastic angular shape",
    }
    records, arrays = [], {}
    with TemporaryDirectory(prefix="li6-separate-") as temporary:
        model, ids = make_model(args.particles, args.batches, temporary)
        try:
            model.init_lib(directory=temporary, output=False, threads=1)
            data = SymmetricTrialData(args.fraction)
            labels["emission_energy"] = (
                f"Li6 MT{data.product['mt']} neutron energy shape"
            )
            selection = dict(
                energy_mt=data.product["mt"],
                energy_product=data.product["product_index"],
                energy_distribution=data.descriptor["distribution_index"],
                xs_mt=2,
                angle_mt=2,
            )
            # A common high-statistics anchor sets the absolute rate scale.
            # Paired differences exclude this shared anchor error; report both.
            openmc.lib.settings.particles = 10 * args.particles
            anchor = run_trial(data, ids, np.zeros(5), 80000)
            openmc.lib.settings.particles = args.particles
            baseline = float(anchor["tritium_rate"][0][0])
            baseline_se = float(anchor["tritium_rate"][1][0])
            arrays["baseline_mesh"] = anchor["tritium_mesh"][0]
            arrays["baseline_depth"] = anchor["tritium_depth"][0]
            for group, indices in groups.items():
                for sign in [-1, 1]:
                    parameters = np.zeros(5)
                    parameters[indices] = sign * args.fraction
                    # Reuse the seed within each perturbed/control pair only.
                    seed = 210000 + len(records)
                    result, details = paired_trial(
                        data, ids, parameters, seed, anchor
                    )
                    key = f'{group}_{"minus" if sign < 0 else "plus"}'
                    for tally in ["tritium_mesh", "tritium_depth"]:
                        arrays[key + "_" + tally] = result[tally][0]
                        arrays[key + "_" + tally + "_paired_se"] = result[
                            tally
                        ][1]
                    records.append(
                        dict(
                            group=group,
                            sign=sign,
                            parameters=parameters.tolist(),
                            seed=seed,
                            **details,
                        )
                    )
                    print(key, details, flush=True)
            data.assert_restored()
        finally:
            model.finalize_lib()
    np.savez_compressed(args.output / "native_responses.npz", **arrays)
    metadata = dict(
        selection=selection,
        baseline=baseline,
        baseline_standard_error=baseline_se,
        fraction=args.fraction,
        particles_per_batch=args.particles,
        batches=args.batches,
        records=records,
        interpretation="Paired errors exclude shared baseline error. Spatial changes are finite-statistics native estimates; no significance mask applied.",
    )
    (args.output / "results.json").write_text(json.dumps(metadata, indent=2))
    middle = SPATIAL_MESH_DIMENSIONS[1] // 2
    base = arrays["baseline_mesh"].reshape(SPATIAL_MESH_DIMENSIONS, order="F")[
        :, middle, :
    ]
    for group in groups:
        changes = []
        for suffix in ["minus", "plus"]:
            mesh = arrays[f"{group}_{suffix}_tritium_mesh"].reshape(
                SPATIAL_MESH_DIMENSIONS, order="F"
            )
            changes.append(mesh[:, middle, :] - base)
        limit = max(np.max(abs(np.asarray(changes))), 1e-30)
        fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
        for ax, values, label in zip(
            axes,
            [base, *changes],
            [
                "Baseline",
                f"−{100*args.fraction:g}% minus baseline",
                f"+{100*args.fraction:g}% minus baseline",
            ],
        ):
            options = (
                {}
                if label == "Baseline"
                else dict(
                    cmap="RdBu_r",
                    norm=TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit),
                )
            )
            pic = ax.pcolormesh(
                np.linspace(-20, 20, 22),
                np.linspace(-20, 20, 22),
                values.T,
                shading="flat",
                **options,
            )
            ax.set(
                title=label,
                xlabel="x [cm]",
                ylabel="z [cm]; beam along +z",
                aspect="equal",
            )
            fig.colorbar(pic, ax=ax, label="(n,t) / voxel / source")
        fig.suptitle(labels[group] + ": central y=0 voxel slice")
        fig.tight_layout()
        fig.savefig(args.output / f"{group}_spatial.png", dpi=180)
        plt.close(fig)
        fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
        depth = np.linspace(1, 39, len(arrays["baseline_depth"]))
        axes[0].plot(
            depth, arrays["baseline_depth"], color="black", label="Baseline"
        )
        subset = [r for r in records if r["group"] == group]
        for r, suffix in zip(subset, ["minus", "plus"]):
            key = f"{group}_{suffix}_tritium_depth"
            axes[0].errorbar(
                depth,
                arrays[key],
                yerr=arrays[key + "_paired_se"],
                label=f'{r["sign"]*100*args.fraction:+g}% shape/scale',
                alpha=0.8,
            )
        axes[0].set(
            xlabel="Depth from beam entrance [cm]",
            ylabel="(n,t) / depth bin / source",
        )
        axes[0].legend()
        # Delta plot cancels common baseline error and makes small effects visible.
        axes[1].errorbar(
            [-100 * args.fraction, 0, 100 * args.fraction],
            [subset[0]["paired_delta"], 0, subset[1]["paired_delta"]],
            yerr=[
                subset[0]["paired_delta_standard_error"],
                0,
                subset[1]["paired_delta_standard_error"],
            ],
            fmt="o-",
            capsize=4,
        )
        axes[1].axhline(0, color="grey", lw=0.8)
        axes[1].set(
            xlabel="Fractional coefficient change [%]",
            ylabel="Change in total (n,t) / source; ±1 SE",
        )
        fig.suptitle(
            labels[group] + "; other perturbation classes fixed at baseline"
        )
        fig.tight_layout()
        fig.savefig(args.output / f"{group}_response.png", dpi=180)
        plt.close(fig)
    print("Saved:", args.output, flush=True)


if __name__ == "__main__":
    main()
