#!/usr/bin/env python3
"""Runtime Li6 sensitivity benchmark and small surrogate training pipeline.

All three nuclear-data perturbations modify native OpenMC objects. These are
controlled sensitivity parameters, not an evaluated nuclear-data covariance
model. A 20 cm radius sphere of pure Li6 at 0.46 g/cm3 is irradiated by an
monoenergetic 14 MeV pencil beam entering along +z. The benchmark is not
a breeder blanket. Tallies are normalized by OpenMC per source neutron.

Example:
    OMP_NUM_THREADS=2 MPLBACKEND=Agg python -m sandbox.examples.fusion_li6_runtime_uq
"""

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory
from copy import deepcopy

from matplotlib.colors import LogNorm, TwoSlopeNorm
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import qmc

import openmc
import openmc.lib

# These five scalar responses become columns in both output CSV files.
# Two additional tallies below retain energy/direction bins for plotting.
# Odd dimensions put the beam axis in the central voxel and y=0 slice.
PERTURBATION_ENERGY_RANGE_EV = (10.0e6, 14.0e6)
XS_MT = 2
# Positive coefficients: larger elastic XS, harder neutron energies and a
# forward angular tilt. Normalization prevents PDFs increasing at every point.
DEMONSTRATION_PARAMETERS = (np.log(3.0), 8.0, 2.0, 3.0, 1.5)
SWEEP_UPPER_BOUNDS = (0.5, 4.0, 2.0, 2.0, 1.0)

SPATIAL_MESH_DIMENSIONS = (21, 21, 21)
SPATIAL_MESH_LOWER = (-20.0, -20.0, -20.0)
SPATIAL_MESH_UPPER = (20.0, 20.0, 20.0)

OBSERVABLES = (
    "tritium_rate",
    "leakage_current",
    "integrated_flux",
    "degraded_leakage",
    "grazing_leakage",
)
# The order is shared by candidate construction, experiment vectors and CSVs.
# a controls XS; b1/b2 control outgoing energy; c1/c2 control scattering angle.
PARAMETERS = (
    "xs_log_scale",
    "energy_linear",
    "energy_quadratic",
    "angle_P1",
    "angle_P2",
)


def energy_basis(x, incident):
    # x contains native outgoing-energy coordinates in eV for one incident table.
    # Scale each table to dimensionless z so the same coefficients are usable
    # across incident energies. The denominator guard also handles zero endpoints.
    # incident is available for energy-dependent bases; this example ignores it.
    z = x / max(float(x.max()), 1.0)
    return np.array([z, z * z])


def angle_basis(mu, incident):
    # P1 tilts the PDF toward forward/backward scattering; P2 changes its shape
    # symmetrically. mu is the native cosine relative to the incoming direction,
    # not an absolute flight direction or a laboratory-frame angle.
    return np.array([mu, (3 * mu * mu - 1) / 2])


def window_mask(energy):
    """Select native incident-energy grid points from 10 through 14 MeV."""
    return (energy >= PERTURBATION_ENERGY_RANGE_EV[0]) & (
        energy <= PERTURBATION_ENERGY_RANGE_EV[1]
    )


def windowed_shape(law, coefficients, basis):
    """Tilt selected incident tables while copying all other tables exactly."""
    result = deepcopy(law)
    selected = np.flatnonzero(window_mask(law["incident_energy"]))
    if not len(selected):
        raise RuntimeError(
            "No native incident tables in the 10–14 MeV interval"
        )
    subset = deepcopy(law)
    subset["incident_energy"] = law["incident_energy"][selected].copy()
    subset["tables"] = [deepcopy(law["tables"][i]) for i in selected]
    changed = openmc.lib.perturb_tabular(subset, coefficients, basis)
    for i, table in zip(selected, changed["tables"]):
        result["tables"][i] = table
    return result


def window_info(grid):
    """Describe control-point window and native interpolation transitions."""
    indices = np.flatnonzero(window_mask(grid))
    if not len(indices):
        raise RuntimeError("No native grid points in requested interval")
    return dict(
        requested_eV=PERTURBATION_ENERGY_RANGE_EV,
        modified_grid_points=len(indices),
        modified_grid_range_eV=[
            float(grid[indices[0]]),
            float(grid[indices[-1]]),
        ],
        adjacent_transition_bounds_eV=[
            float(grid[max(0, indices[0] - 1)]),
            float(grid[min(len(grid) - 1, indices[-1] + 1)]),
        ],
    )


def make_model(particles, batches, directory):
    # A deliberately simple material/geometry isolates nuclear-data responses.
    # The density is a benchmark choice, not a detailed blanket specification.
    material = openmc.Material(name="pure Li6 benchmark")
    material.add_nuclide("Li6", 1.0)
    material.set_density("g/cm3", 0.46)
    # The negative sphere half-space is its interior. Vacuum removes particles
    # on escape; the boundary surface is also our leakage detector.
    outer = openmc.Sphere(r=20.0, boundary_type="vacuum")
    cell = openmc.Cell(fill=material, region=-outer)
    model = openmc.Model(geometry=openmc.Geometry([cell]))
    # Fixed-source transport samples the same external source each batch;
    # there is no eigenvalue iteration or fission-source convergence to manage.
    model.settings.run_mode = "fixed source"
    # A zero-width beam starts just inside the upstream sphere surface. This
    # small inward offset avoids a source exactly on the vacuum boundary;
    # every neutron initially travels along +z across the sphere diameter.
    model.settings.source = openmc.IndependentSource(
        space=openmc.stats.Point((0, 0, -20.0 + 1.0e-6)),
        angle=openmc.stats.Monodirectional((0, 0, 1)),
        energy=openmc.stats.Discrete([14.0e6], [1.0]),
    )
    # particles is histories PER BATCH, so total histories = particles*batches.
    # Batch-to-batch variation supplies the tally standard error of the mean.
    model.settings.particles = particles
    model.settings.batches = batches
    model.settings.verbosity = 1
    # Read tallies through openmc.lib rather than producing a statepoint for
    # every trial. Temporary XML/native output and persistent CSV/PNG outputs
    # have separate directories so repeated runs do not reuse stale files.
    model.settings.output = {
        "path": str(directory),
        "summary": False,
        "tallies": False,
    }
    model.settings.statepoint = {"batches": []}
    neutron = openmc.ParticleFilter(["neutron"])
    surface = openmc.SurfaceFilter([outer])
    cell_filter = openmc.CellFilter([cell])
    # Tally interpretation (all quantities are per source neutron):
    # - (n,t) counts Li6 MT105 reactions; each produces one triton here.
    # - current measures escaped neutron weight, not an escape-event probability;
    #   multiplying reactions can make neutron current exceed one.
    # - cell flux is integrated track length [cm/source], not divided by volume.
    # - E<13 MeV excludes uncollided 14 MeV source neutrons, but includes both
    #   scattered source neutrons and reaction secondaries, not just MT24 daughters.
    # - surface mu is direction dot outward surface normal. Small positive mu
    #   means grazing escape and gives sensitivity to angular-law changes.
    # Surface incidence cosine measures grazing escape independently of the
    # fixed +z direction of the incident pencil beam.
    # Separate (n,t) tallies explain changes in the total production rate:
    # energy bins identify which INCIDENT neutron energies induce reactions;
    # mesh bins show where reactions occur along the beam penetration direction.
    # Mesh scores integrate over each slab's intersection with the sphere.
    tritium_energy_edges = np.array(
        [0.0, 1.0, 1.0e3, 1.0e5, 1.0e6, 5.0e6, 10.0e6, 14.01e6]
    )
    depth_mesh = openmc.RegularMesh()
    depth_mesh.dimension = (1, 1, 20)
    depth_mesh.lower_left = (-20.0, -20.0, -20.0)
    depth_mesh.upper_right = (20.0, 20.0, 20.0)
    # A full 3D mesh resolves lateral spreading as well as beam penetration.
    # The cell filter clips scoring to Li6 inside the sphere. Values are rates
    # integrated in each voxel, including partially filled boundary voxels.
    spatial_mesh = openmc.RegularMesh()
    spatial_mesh.dimension = SPATIAL_MESH_DIMENSIONS
    spatial_mesh.lower_left = SPATIAL_MESH_LOWER
    spatial_mesh.upper_right = SPATIAL_MESH_UPPER
    definitions = [
        ("tritium_rate", [cell_filter, neutron], "(n,t)"),
        (
            "tritium_energy",
            [cell_filter, neutron, openmc.EnergyFilter(tritium_energy_edges)],
            "(n,t)",
        ),
        (
            "tritium_depth",
            [cell_filter, neutron, openmc.MeshFilter(depth_mesh)],
            "(n,t)",
        ),
        (
            "tritium_mesh",
            [cell_filter, neutron, openmc.MeshFilter(spatial_mesh)],
            "(n,t)",
        ),
        ("leakage_current", [surface, neutron], "current"),
        ("integrated_flux", [cell_filter, neutron], "flux"),
        (
            "degraded_leakage",
            [surface, neutron, openmc.EnergyFilter([0.0, 13.0e6])],
            "current",
        ),
        (
            "grazing_leakage",
            [surface, neutron, openmc.MuSurfaceFilter([0.0, 0.5])],
            "current",
        ),
        (
            "leakage_spectrum",
            [
                surface,
                neutron,
                openmc.EnergyFilter(np.linspace(0, 14.01e6, 29)),
            ],
            "current",
        ),
        (
            "leakage_direction",
            [surface, neutron, openmc.MuSurfaceFilter(np.linspace(0, 1, 11))],
            "current",
        ),
    ]
    ids = {}
    for name, filters, score in definitions:
        tally = openmc.Tally(name=name)
        tally.filters = filters
        tally.scores = [score]
        if score == "(n,t)":
            tally.nuclides = ["Li6"]
        # Track-length scores reduce noise for reaction rate and flux; surface
        # crossing current uses the analog estimator. No extra manual source
        # normalization is applied because OpenMC already performs it.
        tally.estimator = "tracklength" if score != "current" else "analog"
        model.tallies.append(tally)
        ids[name] = tally.id
    return model, ids


def equal(a, b):
    # Compare the complete nested runtime snapshots exactly, including PDFs
    # AND loaded CDFs. A tolerance could conceal gradual restoration drift.
    if isinstance(a, dict):
        return a.keys() == b.keys() and all(equal(a[k], b[k]) for k in a)
    if isinstance(a, list):
        return len(a) == len(b) and all(equal(x, y) for x, y in zip(a, b))
    if isinstance(a, np.ndarray):
        return np.array_equal(a, b)
    return a == b


class TrialData:
    """Baseline snapshots and constrained candidates for a loaded Li6 model."""

    def __init__(self):
        # model.init_lib() must already have loaded Li6. These getters copy
        # native C++ data through the C API; no evaluated HDF5 is read in Python.
        self.index = openmc.lib.nuclides["Li6"]._index
        mts = openmc.lib.get_nuclide_mt_numbers(self.index)
        if XS_MT not in mts:
            raise RuntimeError(f"Li6 MT={XS_MT} is absent")
        # Target elastic MT2 XS. The (n,t) tallies remain MT105 observables,
        # so their changes reflect altered transport rather than a direct
        # scaling of the tritium-production reaction cross section.
        self.xs_info = openmc.lib.get_reaction_info("Li6", XS_MT)
        if self.xs_info["redundant"]:
            raise RuntimeError(
                "Cannot perturb a redundant reaction cross section"
            )
        self.grid = openmc.lib.nuclide_energy_grid(self.index)
        self.threshold = openmc.lib.get_reaction_threshold_energy(
            self.index, XS_MT
        )
        self.xs = openmc.lib.get_nuclide_xs(self.index, XS_MT)
        # The XS getter returns only the native array starting at its threshold.
        # The mutation API takes the FULL nuclide grid, so pad below threshold
        # with zeros. Retain the original threshold throughout every trial.
        self.full_xs = np.zeros_like(self.grid)
        self.full_xs[self.grid >= self.threshold] = self.xs
        # Discover the secondary-energy target instead of assuming an MT.
        # Only neutron products with a supported independent energy adapter
        # and an incident grid covering 14 MeV are eligible. Product and
        # distribution indices are zero-based native-vector indices.
        self.products = openmc.lib.get_reaction_products("Li6")
        candidates = []
        for p in self.products:
            if p["particle"] != "neutron":
                continue
            for d in p["distributions"]:
                if d["energy_mutation_supported"]:
                    law = openmc.lib.get_emission_distribution(
                        "Li6",
                        p["mt"],
                        p["product_index"],
                        d["distribution_index"],
                    )["energy"]
                    if (
                        law["incident_energy"][0]
                        <= 14.0e6
                        <= law["incident_energy"][-1]
                    ):
                        candidates.append((p, d, law))
        if not candidates:
            raise RuntimeError("No mutable neutron energy law covering 14 MeV")
        # Use the first eligible native law. In the validated library this is
        # MT24; a different library may expose a different eligible target.
        self.product, self.descriptor, self.energy = candidates[0]
        # Elastic MT2 is the angular demonstration. An isotropic tag has no
        # shape parameter, so this example refuses to convert it into a table.
        self.angle = openmc.lib.get_angular_distribution("Li6", 2)
        if self.angle["type"] != "tabular-angle":
            raise RuntimeError("Elastic angle is not an existing tabular law")
        elastic = openmc.lib.get_emission_distribution("Li6", 2)
        if elastic["type"] != "uncorrelated" or elastic["energy"] is not None:
            raise RuntimeError(
                "Expected angular-only elastic MT2 representation"
            )
        print(
            "MT2 outgoing energy is coupled to angle by native elastic kinematics; "
            "the separate energy-only case uses the discovered nonelastic law."
        )
        for name, grid in (
            ("XS", self.grid),
            ("energy", self.energy["incident_energy"]),
            ("angle", self.angle["incident_energy"]),
        ):
            print(f"{name} 10–14 MeV window: {window_info(grid)}")
        print(
            "Native interpolation transitions extend to the adjacent unchanged "
            "grid points; this is not a hard transport-energy cutoff."
        )
        # REACTION_NAME is only a static MT label lookup; distribution values
        # and native Q/frame/redundancy come from the runtime APIs above.
        print(
            f"XS: Li6 MT={XS_MT} {openmc.data.REACTION_NAME[XS_MT]}, "
            f'native Q={self.xs_info["q_value"]:.9g} eV; '
            f"native threshold grid starts at {self.threshold:.9g} eV; "
            "perturbation control points: 10–14 MeV."
        )
        print(
            f'Energy: Li6 MT={self.product["mt"]} '
            f'{openmc.data.REACTION_NAME[self.product["mt"]]}, '
            f'product={self.product["product_index"]}, '
            f'distribution={self.descriptor["distribution_index"]}, '
            f'frame={self.product["reference_frame"]}; ContinuousTabular'
        )
        print(
            "Angle: Li6 MT=2 product=0 distribution=0; tabular-angle, center-of-mass"
        )

    def candidates(self, parameters):
        # Build independent Python candidates from the SAME stored baseline
        # every time. Preparing these dictionaries does not mutate transport.
        # Omit zero-coefficient operations to retain original native tables
        # exactly, including any small PDF/CDF differences in evaluated data.
        a, b1, b2, c1, c2 = parameters
        args = dict(xs=[], energy=[], angle=[])
        # Constant basis phi(E)=1 gives sigma_new(E)=sigma_original(E)*exp(a).
        # The window indicator is the XS basis: a=log(3) triples MT2
        # only at native incident grid points from 10 through 14 MeV.
        if a:
            args["xs"] = [
                dict(
                    nuclide="Li6",
                    mt=XS_MT,
                    values=openmc.lib.perturb_xs(
                        self.full_xs,
                        [a],
                        window_mask(self.grid).astype(float)[None, :],
                    ),
                )
            ]
        # Tilt the continuous energy PDF by exp(b1*z+b2*z^2). The helper
        # normalizes the resulting native histogram/linear interpolant and
        # preserves grids, interpolation, and any discrete line probabilities.
        if b1 or b2:
            args["energy"] = [
                dict(
                    nuclide="Li6",
                    mt=self.product["mt"],
                    product=self.product["product_index"],
                    distribution=self.descriptor["distribution_index"],
                    law=windowed_shape(self.energy, [b1, b2], energy_basis),
                )
            ]
        # Tilt angular PDFs by exp(c1*P1+c2*P2), normalizing each incident table
        # on -1<=mu<=1. Native interpolation and center-of-mass frame remain.
        if c1 or c2:
            args["angle"] = [
                dict(
                    nuclide="Li6",
                    mt=2,
                    law=windowed_shape(self.angle, [c1, c2], angle_basis),
                )
            ]
        return args

    def assert_restored(self):
        # These are fresh native readbacks, not checks of our cached Python
        # dictionaries. Exact equality verifies restoration without numerical
        # reconstruction of the baseline PDF/CDF arrays.
        assert np.array_equal(
            openmc.lib.get_nuclide_xs(self.index, XS_MT), self.xs
        )
        assert equal(openmc.lib.get_angular_distribution("Li6", 2), self.angle)
        law = openmc.lib.get_emission_distribution(
            "Li6",
            self.product["mt"],
            self.product["product_index"],
            self.descriptor["distribution_index"],
        )["energy"]
        assert equal(law, self.energy)


def run_trial(data, ids, parameters, seed):
    # Detect any persistence from a previous trial before resetting run state.
    data.assert_restored()
    # hard_reset() resets simulation/tally/RNG state; it does NOT restore
    # nuclear data. Restoration is the separate perturbation-context operation.
    openmc.lib.hard_reset()  # clears total_gen, tallies, timers and RNG stride
    openmc.lib.settings.seed = seed
    # The context validates candidates, saves owned native rollback copies,
    # initializes simulation state as needed, applies changes and rebuilds XS.
    # It restores originals on normal exit AND when an exception escapes.
    with openmc.lib.nuclear_data_perturbation(**data.candidates(parameters)):
        # run() uses the actual C++ transport kernel with the modified objects.
        # It finalizes simulation state afterwards but keeps the library/data
        # loaded, allowing the next trial in this same session.
        openmc.lib.run()
        result = {}
        for name, tally_id in ids.items():
            tally = openmc.lib.tallies[tally_id]
            # Copy results now: the next simulation_init() resets accumulation.
            # OpenMC std_dev is the standard error of the tally mean, not the
            # spread of individual particle scores.
            result[name] = (
                tally.mean.ravel().copy(),
                tally.std_dev.ravel().copy(),
            )
            assert tally.num_realizations == openmc.lib.settings.get_batches()
    data.assert_restored()
    return result


def save_figure(fig, path):
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_data(data, directory, parameters):
    # Read the APPLIED native candidates back for plotting. This verifies that
    # plots represent what the transport sampler would use, not just a proposed
    # Python transformation. No transport run is needed to inspect these laws.
    args = data.candidates(parameters)
    with openmc.lib.nuclear_data_perturbation(**args):
        xs = openmc.lib.get_nuclide_xs(data.index, XS_MT)
        energy = openmc.lib.get_emission_distribution(
            "Li6",
            data.product["mt"],
            data.product["product_index"],
            data.descriptor["distribution_index"],
        )["energy"]
        angle = openmc.lib.get_angular_distribution("Li6", 2)
    # Snapshot arrays remain usable after the context restores native data.
    fig, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
    e = data.grid[data.grid >= data.threshold] / 1.0e6
    for values, label in ((data.xs, "Baseline"), (xs, "Perturbed")):
        axes[0].loglog(e, values, label=label)
    # Ratios are undefined at zero baseline XS; leave those points out.
    valid = data.xs > 0
    axes[1].semilogx(e[valid], xs[valid] / data.xs[valid])
    axes[0].set(
        ylabel="Microscopic cross section [barn]",
        title=f"Li6 MT={XS_MT} elastic XS: 10–14 MeV grid window",
    )
    axes[1].set(xlabel="Incident energy [MeV]", ylabel="Perturbed / baseline")
    axes[0].legend()
    # Zoom around the requested window; the full native XS arrays remain stored.
    for ax in axes:
        ax.set_xlim(8.0, 16.0)
        ax.axvspan(10.0, 14.0, color="grey", alpha=0.12)
    save_figure(fig, directory / "cross_section.png")
    for kind, original, changed, key in (
        ("energy", data.energy, energy, "energy"),
        ("angle_mu", data.angle, angle, "x"),
    ):
        fig, ax = plt.subplots(figsize=(8, 5))
        for law, label in ((original, "Baseline"), (changed, "Perturbed")):
            # Plot the nearest actual incident table, reporting its energy.
            # Transport still uses the native incident-energy sampling scheme;
            # this figure does not substitute a nearest-table transport model.
            i = int(np.argmin(abs(law["incident_energy"] - 14.0e6)))
            table = law["tables"][i]
            # The first n energy entries are discrete probability masses.
            # Angular Tabular objects have no such atoms, hence the default 0.
            n = table.get("n_discrete", 0)
            # For energy: x[eV]/1e6 -> MeV, p[1/eV]*1e6 -> 1/MeV.
            # mu is dimensionless and its PDF needs no degree/radian Jacobian.
            scale = 1.0e6 if key == "energy" else 1.0
            ax.plot(
                table[key][n:] / scale,
                table["pdf"][n:] * scale,
                label=label,
                drawstyle=(
                    "steps-post" if table["interpolation"] == 1 else "default"
                ),
            )
            # Atoms are probabilities, not densities: do not multiply by an
            # energy Jacobian or draw them as part of the continuous PDF.
            if n:
                atoms = ax.twinx()
                atoms.vlines(table[key][:n] / scale, 0, table["pdf"][:n])
                atoms.set_ylabel("Discrete probability mass")
            print(
                f"{kind} {label}: selected incident energy "
                f'{law["incident_energy"][i]:.9g} eV; interpolation '
                f'{table["interpolation"]}; discrete count {n}'
            )
        ax.legend()
        if key == "energy":
            ax.set(
                xlabel="Outgoing neutron energy [MeV]",
                ylabel="PDF [1/MeV]",
                title=f'Li6 MT={data.product["mt"]}, {data.product["reference_frame"]}, 14 MeV table',
            )
        else:
            ax.set(
                xlabel="Scattering cosine, μ = cos(θ)",
                ylabel="Probability density in μ",
                xlim=(-1, 1),
                title="Li6 MT=2, center-of-mass, 14 MeV table",
            )
        save_figure(fig, directory / f"{kind}.png")


def write_rows(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def row(name, parameters, seed, particles, batches, result=None, error=""):
    # One CSV row contains inputs, run settings and scalar responses/errors.
    # Failed samples keep their inputs and an error message, with blank targets
    # rather than fabricated zero-valued observations for surrogate training.
    values = dict(
        case=name,
        **dict(zip(PARAMETERS, map(float, parameters))),
        seed=seed,
        particles_per_batch=particles,
        batches=batches,
        histories=particles * batches,
        nuclide="Li6",
        xs_mt=XS_MT,
        angle_mt=2,
        angle_product=0,
        angle_distribution=0,
        status="failed" if error else "success",
        error=error,
    )
    for observable in OBSERVABLES:
        values[observable + "_mean"] = (
            float(result[observable][0][0]) if result else ""
        )
        values[observable + "_std_dev"] = (
            float(result[observable][1][0]) if result else ""
        )
    return values


def plot_tritium_comparisons(results, cases, directory):
    """Save total, energy-resolved and depth-resolved (n,t) comparisons."""
    baseline, baseline_sd = results["Baseline"]["tritium_rate"]
    comparisons, bin_rows = [], []
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for i, (name, _, _) in enumerate(cases):
        mean, sd = results[name]["tritium_rate"]
        # Independent case seeds permit quadrature combination of standard
        # errors. Baseline/restored baseline share histories and are compared
        # exactly, so their difference and its uncertainty are zero.
        same_histories = name in ("Baseline", "Restored baseline")
        change = 100.0 * (mean[0] / baseline[0] - 1.0)
        change_sd = (
            0.0
            if same_histories
            else 100.0
            * np.hypot(
                sd[0] / baseline[0],
                mean[0] * baseline_sd[0] / baseline[0] ** 2,
            )
        )
        axes[0].errorbar(i, mean[0], yerr=sd[0], fmt="o", capsize=4)
        axes[1].errorbar(i, change, yerr=change_sd, fmt="o", capsize=4)
        comparisons.append(
            dict(
                case=name,
                rate_per_source=mean[0],
                std_error=sd[0],
                change_percent=change,
                change_std_error_percent=change_sd,
            )
        )
        print(
            f"(n,t) {name}: {mean[0]:.6g} +/- {sd[0]:.3g} per source; "
            f"change = {change:+.2f}% +/- {change_sd:.2f}%"
        )
    for ax in axes:
        ax.set_xticks(
            range(len(cases)), [c[0] for c in cases], rotation=45, ha="right"
        )
        ax.grid(axis="y", alpha=0.25)
    axes[0].set_ylabel("Li6 (n,t) reactions / source neutron")
    axes[1].set_ylabel("Change from baseline [%]")
    axes[1].axhline(0.0, color="grey", linewidth=0.8)
    fig.suptitle(
        "14 MeV pencil beam: Li6 tritium production (one standard error)"
    )
    save_figure(fig, directory / "tritium_comparison.png")
    write_rows(directory / "tritium_comparison.csv", comparisons)

    energy_edges = np.array(
        [0.0, 1.0, 1.0e3, 1.0e5, 1.0e6, 5.0e6, 10.0e6, 14.01e6]
    )
    depth_edges = np.linspace(-20.0, 20.0, 21)
    for tally_name, edges, xlabel in (
        ("tritium_energy", energy_edges, "Incident neutron energy [eV]"),
        ("tritium_depth", depth_edges, "Depth from beam entrance [cm]"),
    ):
        fig, ax = plt.subplots(figsize=(10, 5))
        centres = (edges[:-1] + edges[1:]) / 2
        if tally_name == "tritium_depth":
            centres += 20.0
        else:
            centres[1:] = np.sqrt(edges[1:-1] * edges[2:])
            ax.set_xscale("log")
        for name, _, _ in cases:
            mean, sd = results[name][tally_name]
            ax.errorbar(centres, mean, yerr=sd, label=name, capsize=2)
            for j, (value, error) in enumerate(zip(mean, sd)):
                bin_rows.append(
                    dict(
                        case=name,
                        tally=tally_name,
                        lower=edges[j],
                        upper=edges[j + 1],
                        coordinate_units=(
                            "eV"
                            if tally_name == "tritium_energy"
                            else "z [cm]"
                        ),
                        rate_per_source=value,
                        std_error=error,
                    )
                )
        # These are bin-integrated rates, not PDFs or rates per unit energy.
        ax.set(
            xlabel=xlabel, ylabel="Li6 (n,t) reactions / bin / source neutron"
        )
        ax.legend(fontsize=8)
        ax.grid(alpha=0.2)
        save_figure(fig, directory / f"{tally_name}.png")
    write_rows(directory / "tritium_bins.csv", bin_rows)


def plot_tritium_mesh(results, cases, directory):
    """Save full 3D rates/errors and central x-z slice comparisons."""
    names = [case[0] for case in cases]
    # RegularMesh bins enumerate x fastest, then y, then z. Fortran order
    # maps flattened native bins to arrays indexed [x, y, z] correctly.
    means = np.stack(
        [
            results[name]["tritium_mesh"][0].reshape(
                SPATIAL_MESH_DIMENSIONS, order="F"
            )
            for name in names
        ]
    )
    errors = np.stack(
        [
            results[name]["tritium_mesh"][1].reshape(
                SPATIAL_MESH_DIMENSIONS, order="F"
            )
            for name in names
        ]
    )
    edges = [
        np.linspace(lo, hi, n + 1)
        for lo, hi, n in zip(
            SPATIAL_MESH_LOWER, SPATIAL_MESH_UPPER, SPATIAL_MESH_DIMENSIONS
        )
    ]
    # Keep all voxels and native standard errors for downstream 3D analysis.
    # Zero-score bins may have infinite reported errors in the runtime API;
    # preserve these rather than manufacturing certainty for unsampled regions.
    np.savez_compressed(
        directory / "tritium_mesh.npz",
        cases=np.array(names),
        mean=means,
        std_error=errors,
        x_edges_cm=edges[0],
        y_edges_cm=edges[1],
        z_edges_cm=edges[2],
    )
    middle = SPATIAL_MESH_DIMENSIONS[1] // 2
    slices = means[:, :, middle, :]
    # This is one finite-thickness central voxel slice, not a sum over y.
    # Consequently there is no invalid quadrature sum of correlated bin errors.
    positive = slices[slices > 0]
    norm = LogNorm(vmin=positive.min(), vmax=positive.max())
    differences = slices - slices[0]
    maximum = max(float(np.max(abs(differences))), np.finfo(float).tiny)
    difference_norm = TwoSlopeNorm(vmin=-maximum, vcenter=0.0, vmax=maximum)
    for filename, arrays, color_norm, cmap, label in (
        (
            "tritium_mesh_slices.png",
            slices,
            norm,
            "viridis",
            "Li6 (n,t) reactions / voxel / source neutron",
        ),
        (
            "tritium_mesh_changes.png",
            differences,
            difference_norm,
            "RdBu_r",
            "Change from baseline [reactions / voxel / source neutron]",
        ),
    ):
        fig, axes = plt.subplots(
            2,
            3,
            figsize=(13, 9),
            sharex=True,
            sharey=True,
            layout="constrained",
        )
        for i, ax in enumerate(axes.ravel()):
            values = arrays[i]
            if filename == "tritium_mesh_slices.png":
                values = np.ma.masked_less_equal(values, 0.0)
            # pcolormesh expects rows along z and columns along x.
            picture = ax.pcolormesh(
                edges[0],
                edges[2],
                values.T,
                norm=color_norm,
                cmap=cmap,
                shading="flat",
            )
            ax.set_title(names[i])
            ax.set_aspect("equal")
            ax.set_xlabel("x [cm]")
            ax.set_ylabel("z [cm]; beam travels along +z")
            # Mark the material boundary and beam entrance for orientation.
            ax.add_patch(
                plt.Circle(
                    (0.0, 0.0), 20.0, fill=False, color="grey", linewidth=0.8
                )
            )
            ax.annotate(
                "",
                xy=(0.0, -16.0),
                xytext=(0.0, -20.0),
                arrowprops=dict(arrowstyle="->", color="black"),
            )
        fig.colorbar(
            picture, ax=axes.ravel().tolist(), label=label, shrink=0.85
        )
        fig.suptitle(
            "14 MeV pencil beam: central y=0 mesh slice "
            f"({edges[1][middle+1]-edges[1][middle]:.3f} cm thick)"
        )
        # Constrained layout accommodates the shared colorbar; applying
        # tight_layout here would disable it and distort the panel sizes.
        fig.savefig(directory / filename, dpi=180)
        plt.close(fig)
    print(
        "Saved tritium_mesh_slices.png, tritium_mesh_changes.png and "
        "tritium_mesh.npz (full 3D means and standard errors)."
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--particles", type=int, default=20000)
    parser.add_argument("--batches", type=int, default=20)
    parser.add_argument("--sweep-samples", type=int, default=8)
    parser.add_argument("--sweep-particles", type=int, default=5000)
    parser.add_argument("--stress-cycles", type=int, default=50)
    parser.add_argument(
        "--output", type=Path, default=Path("sandbox/examples/fusion_li6_uq_output")
    )
    args = parser.parse_args()
    if (
        min(args.particles, args.sweep_particles) < 1
        or args.batches < 2
        or args.stress_cycles < 0
        or args.sweep_samples < 0
    ):
        parser.error(
            "positive particles, at least 2 batches, nonnegative stress cycles required"
        )
    args.output.mkdir(parents=True, exist_ok=True)
    # Deliberately strong, admissible sensitivity changes make responses visible
    # with manageable statistics. These are not sampled evaluated covariances.
    combined = list(DEMONSTRATION_PARAMETERS)
    # Perturbed cases use independent seeds for interpretable comparisons.
    # Only baseline/restored baseline share a seed, to test exact reproduction.
    cases = [
        ("Baseline", [0.0] * 5, 91),
        ("XS only", [combined[0], 0, 0, 0, 0], 92),
        ("Energy only", [0, *combined[1:3], 0, 0], 93),
        ("Angle only", [0, 0, 0, *combined[3:]], 94),
        ("All three", combined, 95),
        ("Restored baseline", [0.0] * 5, 91),
    ]
    # Initialize the runtime library ONCE for controlled cases, the sweep and
    # stress test. Temporary input XML is removed only after library finalization.
    with TemporaryDirectory(prefix="openmc-li6-uq-") as temporary:
        model, ids = make_model(args.particles, args.batches, temporary)
        try:
            # A single transport thread avoids last-bit variation from parallel
            # tally reductions in the deterministic restoration comparisons.
            model.init_lib(directory=temporary, output=False, threads=1)
            data = TrialData()
            plot_data(data, args.output, combined)
            # Plot inspection entered a mutation context and initialized
            # simulation state. Close that state before starting the first run.
            openmc.lib.simulation_finalize()
            results, rows = {}, []
            for name, parameters, seed in cases:
                result = run_trial(data, ids, parameters, seed)
                results[name] = result
                record = row(
                    name,
                    parameters,
                    seed,
                    args.particles,
                    args.batches,
                    result,
                )
                record.update(
                    energy_mt=data.product["mt"],
                    energy_product=data.product["product_index"],
                    energy_distribution=data.descriptor["distribution_index"],
                )
                rows.append(record)
                print(
                    name,
                    {
                        o: (record[o + "_mean"], record[o + "_std_dev"])
                        for o in OBSERVABLES
                    },
                    flush=True,
                )
            # Check all tally means AND standard errors, including spectrum
            # and directional bins, against the same-seed restored baseline.
            for name in ids:
                for a, b in zip(
                    results["Baseline"][name],
                    results["Restored baseline"][name],
                ):
                    np.testing.assert_array_equal(a, b)
            print(
                "Restored baseline tallies and uncertainties are bitwise identical.",
                flush=True,
            )
            write_rows(args.output / "controlled_cases.csv", rows)
            plot_tritium_comparisons(results, cases, args.output)
            plot_tritium_mesh(results, cases, args.output)
            # Error bars are one standard error. Small effects must be assessed
            # against their uncertainties, not dismissed from similar means.
            fig, axes = plt.subplots(2, 3, figsize=(14, 8))
            for observable, ax in zip(OBSERVABLES, axes.ravel()):
                means = [results[n][observable][0][0] for n, *_ in cases]
                sd = [results[n][observable][1][0] for n, *_ in cases]
                ax.errorbar(np.arange(6), means, yerr=sd, fmt="o", capsize=3)
                ax.set_xticks(
                    np.arange(6),
                    [n for n, *_ in cases],
                    rotation=45,
                    ha="right",
                )
                ax.set_title(observable.replace("_", " "))
                ax.set_ylabel(
                    "Per source neutron"
                    if observable != "integrated_flux"
                    else "Track length [cm/source neutron]"
                )
            axes.ravel()[-1].axis("off")
            save_figure(fig, args.output / "transport_responses.png")
            for name, edges, label in (
                (
                    "leakage_spectrum",
                    np.linspace(0, 14.01, 29),
                    "Escaping neutron energy [MeV]",
                ),
                (
                    "leakage_direction",
                    np.linspace(0, 1, 11),
                    "Outward incidence cosine at sphere surface",
                ),
            ):
                fig, ax = plt.subplots(figsize=(9, 5))
                centres = (edges[:-1] + edges[1:]) / 2
                for case, result in results.items():
                    ax.errorbar(
                        centres,
                        result[name][0],
                        yerr=result[name][1],
                        label=case,
                        capsize=2,
                    )
                ax.set(
                    xlabel=label,
                    ylabel="Escaped neutron current per bin/source neutron",
                )
                ax.legend(fontsize=8)
                save_figure(fig, args.output / f"{name}.png")
            # Stage 6: a deterministic five-dimensional Latin hypercube.
            design = qmc.LatinHypercube(5, seed=20261009).random(
                args.sweep_samples
            )
            # Convert unit-cube coordinates to the declared coefficient ranges:
            # All coefficients are nonnegative: a in [0,.5], b1 in [0,4],
            # b2 in [0,2], c1 in [0,2], c2 in [0,1]. This is an experimental design,
            # not an independent nuclear-data covariance sampling model.
            design = design * np.array(SWEEP_UPPER_BOUNDS)
            # Use fewer histories for this demonstration training dataset.
            # run_trial still resets state and restores native data for EVERY row.
            openmc.lib.settings.particles = args.sweep_particles
            training = []
            for i, parameters in enumerate(design):
                seed = 1000 + i
                error, result = "", None
                try:
                    result = run_trial(data, ids, parameters, seed)
                except Exception as problem:
                    # Record candidate/run failure, but do not continue if the
                    # original nuclear data could not be restored: the assertion
                    # below then raises and stops generation of invalid samples.
                    error = str(problem)
                    data.assert_restored()
                record = row(
                    f"sample_{i}",
                    parameters,
                    seed,
                    args.sweep_particles,
                    args.batches,
                    result,
                    error,
                )
                record.update(
                    energy_mt=data.product["mt"],
                    energy_product=data.product["product_index"],
                    energy_distribution=data.descriptor["distribution_index"],
                )
                training.append(record)
            if training:
                write_rows(args.output / "training.csv", training)
            # Stage 7 reliability: alternating combined trials and exact baseline.
            openmc.lib.settings.particles = 200
            openmc.lib.settings.set_batches(2, add_sp_batch=False)
            # Small histories make many restoration cycles affordable. This
            # stress test checks determinism/state safety, not response precision.
            # Every cycle includes a perturbed run and a same-seed baseline run.
            reference = run_trial(data, ids, [0.0] * 5, 713)
            rss = []
            for i in range(args.stress_cycles):
                run_trial(data, ids, combined, 712)
                restored = run_trial(data, ids, [0.0] * 5, 713)
                for name in ids:
                    for a, b in zip(reference[name], restored[name]):
                        np.testing.assert_array_equal(a, b)
                # Linux resident pages give a coarse memory-growth diagnostic.
                # A plateau is useful evidence; it is not a full leak-detector test.
                rss.append(
                    int(Path("/proc/self/statm").read_text().split()[1])
                    * os.sysconf("SC_PAGE_SIZE")
                )
            stress = dict(
                cycles=args.stress_cycles,
                transport_runs=2 * args.stress_cycles,
                exact_restoration=True,
                exact_baseline_reproducibility=True,
                resident_bytes=rss,
                rss_growth_bytes=rss[-1] - rss[0] if rss else 0,
            )
            (args.output / "stress.json").write_text(
                json.dumps(stress, indent=2)
            )
            # Read only the cross_sections.xml INDEX for its identifier/hash.
            # No evaluated HDF5 distribution data are obtained through this path.
            library = Path(os.environ["OPENMC_CROSS_SECTIONS"])
            # Metadata supplies the basis/physics conventions and provenance
            # needed to interpret and reproduce the tabular training rows.
            metadata = dict(
                openmc_version=openmc.__version__,
                git_commit=subprocess.check_output(
                    ["git", "rev-parse", "HEAD"], text=True
                ).strip(),
                git_dirty=True,
                transport_threads=1,
                library=str(library),
                library_index_sha256=hashlib.sha256(
                    library.read_bytes()
                ).hexdigest(),
                tritium_mesh=dict(
                    dimensions=SPATIAL_MESH_DIMENSIONS,
                    lower_left_cm=SPATIAL_MESH_LOWER,
                    upper_right_cm=SPATIAL_MESH_UPPER,
                    plot="Central y=0 voxel slice, finite thickness; no projection",
                    units="reactions per voxel per source neutron",
                ),
                geometry=dict(
                    radius_cm=20.0, density_g_cm3=0.46, material="pure Li6"
                ),
                source=dict(
                    energy_eV=14.0e6,
                    position=[0, 0, -20.0 + 1.0e-6],
                    angle="monodirectional",
                    direction=[0, 0, 1],
                    description="Pencil beam just inside upstream surface",
                ),
                selection=dict(
                    xs_mt=XS_MT,
                    xs_metadata=data.xs_info,
                    energy_mt=data.product["mt"],
                    energy_product=data.product["product_index"],
                    energy_distribution=data.descriptor["distribution_index"],
                    angle_mt=2,
                ),
                incident_energy_window=dict(
                    xs=window_info(data.grid),
                    energy=window_info(data.energy["incident_energy"]),
                    angle=window_info(data.angle["incident_energy"]),
                    semantics="Only in-window native control points/tables change. Native interpolation creates transitions to adjacent unchanged points outside the window.",
                    elastic_energy="MT2 has no independent energy PDF; angular changes alter outgoing energies through native elastic kinematics.",
                ),
                basis=dict(
                    xs="phi(E)=1 at native grid points 10–14 MeV, else 0; sigma*=exp(a*phi)",
                    energy="z=E_out/max(native table E_out); psi=[z,z^2]",
                    angle="P1=mu, P2=(3*mu^2-1)/2",
                ),
                shape="Exponential tilt at native points, native interpolant normalized; discrete masses retained.",
                controlled_parameters=dict(zip(PARAMETERS, combined)),
                perturbation_direction="Increased MT2 XS, harder MT24 outgoing energy and forward MT2 angular tilt. PDFs remain normalized; response tallies may increase or decrease.",
                design=dict(
                    lower_bounds=[0.0] * 5,
                    upper_bounds=SWEEP_UPPER_BOUNDS,
                    type="Latin hypercube",
                    seed=20261009,
                    samples=args.sweep_samples,
                ),
                comparison="Independent seeds for perturbation cases; seed 91 shared only for exact restored baseline.",
                observables=dict(
                    tritium_rate="Li6 MT105 reaction rate/source; one triton per (n,t), not a general yield tally",
                    tritium_energy="MT105 rate/source binned by incident neutron energy, not triton energy",
                    tritium_depth="MT105 rate/source integrated in 2 cm slabs along +z",
                    tritium_mesh="MT105 rate/source integrated per 3D voxel; x fastest (Fortran reshape)",
                    leakage_current="outward neutron current/source, may exceed unity from multiplication",
                    integrated_flux="volume-integrated track length cm/source; not divided by sphere volume",
                    degraded_leakage="outward current with E<13 MeV; includes elastic slowing and reaction secondaries",
                    grazing_leakage="outward current with surface incidence cosine 0..0.5",
                ),
                uncertainty="OpenMC standard error across batches; plots show one standard error.",
                limitations="Sensitivity design, not evaluated covariance; only independent tabular law mutation.",
                stress=stress,
            )
            (args.output / "metadata.json").write_text(
                json.dumps(metadata, indent=2)
            )
            print(
                "Stress:",
                {k: v for k, v in stress.items() if k != "resident_bytes"},
            )
            print("Output:", args.output)
        finally:
            # Release the library once, even if a run, plot or validation fails.
            # Never finalize it INSIDE an active perturbation context: native
            # rollback handles refer to objects owned by the loaded session.
            model.finalize_lib()


if __name__ == "__main__":
    main()
