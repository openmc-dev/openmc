"""Zero-centering, normalization and +/-10% bounds of the new scan prior."""

from itertools import product
from copy import deepcopy
import numpy as np
import pytest

from sandbox.examples.li6_symmetric_perturbations import symmetric_shape, paired_design
from sandbox.examples.fusion_li6_runtime_uq import energy_basis, angle_basis
from openmc.lib.perturbation import _areas


@pytest.mark.parametrize("interpolation", [1, 2])
@pytest.mark.parametrize("kind", ["energy", "angle"])
@pytest.mark.parametrize("limit", [0.1, 0.2])
def test_symmetric_pdf_center_and_bounds(interpolation, kind, limit):
    if kind == "energy":
        table = dict(
            energy=np.array([0.5, 0.0, 1.0, 2.0]),
            pdf=np.array([0.2, 0.4, 0.4, 0.4]),
            cdf=np.array([0.2, 0.2, 0.6, 1.0]),
            interpolation=interpolation,
            n_discrete=1,
        )
        law = dict(
            type="continuous-tabular",
            incident_energy=np.array([9.0e6, 10.0e6, 14.0e6, 15.0e6]),
            tables=[deepcopy(table) for _ in range(4)],
        )
        basis, key, n = energy_basis, "energy", 1
    else:
        p = (
            np.array([0.25, 0.75, 0.75])
            if interpolation == 1
            else np.array([0.25, 0.5, 0.75])
        )
        x = np.array([-1.0, 0.0, 1.0])
        table = dict(
            x=x,
            pdf=p,
            cdf=np.r_[0.0, np.cumsum(_areas(x, p, interpolation))],
            type="tabular",
            interpolation=interpolation,
            integral=1.0,
        )
        law = dict(
            type="tabular-angle",
            incident_energy=np.array([9.0e6, 10.0e6, 14.0e6, 15.0e6]),
            tables=[deepcopy(table) for _ in range(4)],
        )
        basis, key, n = angle_basis, "x", 0
    original = table["pdf"].copy()
    for fractions in product([-limit, limit], repeat=2):
        up = symmetric_shape(law, fractions, basis, limit)
        down = symmetric_shape(law, -np.array(fractions), basis, limit)
        for i in [0, 3]:
            np.testing.assert_array_equal(up["tables"][i]["pdf"], original)
            np.testing.assert_array_equal(
                up["tables"][i]["cdf"], law["tables"][i]["cdf"]
            )
        p = up["tables"][1]["pdf"]
        np.testing.assert_allclose(
            (p + down["tables"][1]["pdf"]) / 2, original, atol=1.0e-15
        )
        assert np.all(p >= (1 - limit) * original - 1.0e-15)
        assert np.all(p <= (1 + limit) * original + 1.0e-15)
        assert p[:n].tolist() == original[:n].tolist()
        assert sum(p[:n]) + sum(
            _areas(table[key][n:], p[n:], interpolation)
        ) == pytest.approx(1.0)
        np.testing.assert_allclose(up["tables"][1]["cdf"][-1], 1.0)
    np.testing.assert_array_equal(law["tables"][1]["pdf"], original)


def test_design_exactly_balanced():
    units = paired_design(128, 59)
    fractions = -0.1 + 0.2 * units
    np.testing.assert_allclose(fractions[:64], -fractions[64:], atol=3.0e-17)
    np.testing.assert_allclose(fractions.mean(axis=0), 0.0, atol=1.0e-17)
    assert np.all(abs(fractions) <= 0.1)
    # Linear XS multiplier is centered at one, unlike exp(symmetric coefficient).
    np.testing.assert_allclose((1.0 + fractions[:, 0]).mean(), 1.0)


@pytest.mark.parametrize("limit", [0.1, 0.2])
def test_native_symmetric_bounds(li6, limit):
    import openmc.lib
    from sandbox.examples.li6_symmetric_perturbations import SymmetricTrialData
    from sandbox.examples.fusion_li6_runtime_uq import equal, window_mask

    data = SymmetricTrialData(limit)
    unchanged_mt105 = openmc.lib.get_nuclide_xs(data.index, 105)
    for sign in [-1, 1]:
        fractions = np.full(5, sign * limit)
        with openmc.lib.nuclear_data_perturbation(
            **data.candidates(fractions)
        ):
            actual = openmc.lib.get_nuclide_xs(data.index, 2)
            mask = window_mask(data.grid[data.grid >= data.threshold])
            np.testing.assert_allclose(
                actual[mask], data.xs[mask] * (1 + sign * limit)
            )
            for before, after, key in (
                (
                    data.angle,
                    openmc.lib.get_angular_distribution("Li6", 2),
                    "x",
                ),
                (
                    data.energy,
                    openmc.lib.get_emission_distribution(
                        "Li6",
                        data.product["mt"],
                        data.product["product_index"],
                        data.descriptor["distribution_index"],
                    )["energy"],
                    "energy",
                ),
            ):
                for i, (a, b) in enumerate(
                    zip(before["tables"], after["tables"])
                ):
                    if not window_mask(before["incident_energy"])[i]:
                        assert equal(a, b)
                        continue
                    positive = a["pdf"] > 0
                    ratios = b["pdf"][positive] / a["pdf"][positive]
                    assert np.all(ratios >= (1 - limit) - 1.0e-6)
                    assert np.all(ratios <= (1 + limit) + 1.0e-6)
                    n = b.get("n_discrete", 0)
                    mass = sum(b["pdf"][:n]) + sum(
                        _areas(b[key][n:], b["pdf"][n:], b["interpolation"])
                    )
                    assert mass == pytest.approx(1.0, abs=1.0e-13)
            np.testing.assert_array_equal(
                openmc.lib.get_nuclide_xs(data.index, 105), unchanged_mt105
            )
        data.assert_restored()


def test_batch_estimates_match_native_tally_means(tmp_path):
    import openmc.lib
    from sandbox.examples.fusion_li6_runtime_uq import make_model, run_trial
    from sandbox.examples.li6_symmetric_perturbations import (
        SymmetricTrialData,
        batch_trial,
    )

    model, ids = make_model(300, 4, str(tmp_path))
    try:
        model.init_lib(directory=str(tmp_path), output=False, threads=1)
        data = SymmetricTrialData()
        expected = run_trial(data, ids, [0.0] * 5, 19)
        batches = batch_trial(data, ids, [0.0] * 5, 19)
        for name in ids:
            np.testing.assert_allclose(
                batches[name].mean(axis=0),
                expected[name][0],
                rtol=1.0e-14,
                atol=1.0e-16,
            )
    finally:
        model.finalize_lib()


from tests.unit_tests.test_runtime_uq import li6
