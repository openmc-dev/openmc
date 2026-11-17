"""Matched native-data restoration plus independent transport response checks."""

import numpy as np

import openmc.lib
from sandbox.examples.fusion_li6_runtime_uq import TrialData, make_model, run_trial


def test_runtime_nuclear_data_uq(tmp_path):
    from openmc.mpi import comm
    from tests.regression_tests import config

    intracomm = comm if config["mpi"] else None
    model, ids = make_model(10000, 10, tmp_path)
    assert isinstance(
        model.settings.source[0].angle, openmc.stats.Monodirectional
    )
    np.testing.assert_array_equal(
        model.settings.source[0].angle.reference_uvw, [0, 0, 1]
    )
    try:
        model.init_lib(
            directory=str(tmp_path),
            output=False,
            threads=1,
            intracomm=intracomm,
        )
        data = TrialData()
        baseline = run_trial(data, ids, [0.0] * 5, 91)
        # Both binning schemes partition the same Li6 (n,t) rate, so their
        # integrated means must agree with the unfiltered tally.
        for name in ("tritium_energy", "tritium_depth", "tritium_mesh"):
            np.testing.assert_allclose(
                baseline[name][0].sum(),
                baseline["tritium_rate"][0][0],
                rtol=1.0e-12,
            )
        parameters = [
            [np.log(3.0), 0, 0, 0, 0],
            [0, -8, 2, 0, 0],
            [0, 0, 0, -3, 1.5],
            [np.log(3.0), -8, 2, -3, 1.5],
        ]
        # Elastic XS does not directly scale MT105 production anymore; use
        # degraded leakage to measure its strong transport effect. With angular
        # changes restricted to fast neutrons, integrated flux is a less noisy
        # response than grazing leakage at this test's particle count.
        observables = [
            "degraded_leakage",
            "degraded_leakage",
            "integrated_flux",
            "tritium_rate",
        ]
        for i, (p, observable) in enumerate(zip(parameters, observables)):
            result = run_trial(data, ids, p, 92 + i)
            diff = abs(result[observable][0][0] - baseline[observable][0][0])
            se = np.hypot(result[observable][1][0], baseline[observable][1][0])
            assert diff > 3 * se
        restored = run_trial(data, ids, [0.0] * 5, 91)
        for name in ids:
            for original, actual in zip(baseline[name], restored[name]):
                np.testing.assert_array_equal(original, actual)
    finally:
        model.finalize_lib()


def test_incident_energy_window(tmp_path):
    from sandbox.examples.fusion_li6_runtime_uq import window_mask, equal, XS_MT

    model, ids = make_model(100, 2, tmp_path)
    try:
        model.init_lib(directory=str(tmp_path), output=False, threads=1)
        data = TrialData()
        assert XS_MT == 2
        arguments = data.candidates([np.log(3.0), -8.0, 2.0, -3.0, 1.5])
        with openmc.lib.nuclear_data_perturbation(**arguments):
            changed = openmc.lib.get_nuclide_xs(data.index, XS_MT)
            grid = data.grid[data.grid >= data.threshold]
            np.testing.assert_array_equal(
                changed[~window_mask(grid)], data.xs[~window_mask(grid)]
            )
            np.testing.assert_allclose(
                changed[window_mask(grid)], data.xs[window_mask(grid)] * 3.0
            )
            for baseline, native in (
                (data.angle, openmc.lib.get_angular_distribution("Li6", 2)),
                (
                    data.energy,
                    openmc.lib.get_emission_distribution(
                        "Li6",
                        data.product["mt"],
                        data.product["product_index"],
                        data.descriptor["distribution_index"],
                    )["energy"],
                ),
            ):
                mask = window_mask(baseline["incident_energy"])
                for i in np.flatnonzero(~mask):
                    assert equal(native["tables"][i], baseline["tables"][i])
                assert any(
                    not equal(native["tables"][i], baseline["tables"][i])
                    for i in np.flatnonzero(mask)
                )
        data.assert_restored()
    finally:
        model.finalize_lib()


def test_default_perturbations_increase(tmp_path):
    from sandbox.examples.fusion_li6_runtime_uq import (
        DEMONSTRATION_PARAMETERS,
        XS_MT,
        window_mask,
    )

    model, ids = make_model(100, 2, tmp_path)
    try:
        model.init_lib(directory=str(tmp_path), output=False, threads=1)
        data = TrialData()
        mt, pi, di = (
            data.product["mt"],
            data.product["product_index"],
            data.descriptor["distribution_index"],
        )
        before_mu = openmc.lib.sample_angular_distribution(
            "Li6", 2, 14.0e6, 30000, seed=59
        )
        before_energy = openmc.lib.sample_emission_energy(
            "Li6", mt, 14.0e6, 30000, pi, di, seed=59
        )
        with openmc.lib.nuclear_data_perturbation(
            **data.candidates(DEMONSTRATION_PARAMETERS)
        ):
            values = openmc.lib.get_nuclide_xs(data.index, XS_MT)
            mask = window_mask(data.grid[data.grid >= data.threshold])
            np.testing.assert_allclose(values[mask], data.xs[mask] * 3.0)
            after_mu = openmc.lib.sample_angular_distribution(
                "Li6", 2, 14.0e6, 30000, seed=59
            )
            after_energy = openmc.lib.sample_emission_energy(
                "Li6", mt, 14.0e6, 30000, pi, di, seed=59
            )
            assert after_mu.mean() > before_mu.mean() + 0.05
            assert after_energy.mean() > before_energy.mean() * 1.2
        data.assert_restored()
    finally:
        model.finalize_lib()
