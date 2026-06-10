"""Tests for the ``weight_windows_file`` argument of
:meth:`openmc.Model.convert_to_multigroup`.

The bootstrapped generation itself is covered by the
``random_ray_auto_convert_bootstrap`` regression test; the tests here only
cover the argument validation, which requires no transport simulation or
nuclear data.
"""

from pathlib import Path

import pytest

import openmc


def _minimal_model():
    openmc.reset_auto_ids()
    model = openmc.Model()

    mat = openmc.Material(name="fuel")
    mat.add_nuclide("U235", 1.0)
    mat.set_density("g/cm3", 1.0)
    model.materials.append(mat)

    sph = openmc.Sphere(r=10.0, boundary_type="vacuum")
    model.geometry = openmc.Geometry([openmc.Cell(fill=mat, region=-sph)])

    model.settings.run_mode = "fixed source"
    model.settings.particles = 100
    model.settings.batches = 10
    model.settings.source = openmc.IndependentSource(space=openmc.stats.Point())
    return model


@pytest.mark.parametrize("method", ["stochastic_slab", "infinite_medium"])
def test_weight_windows_file_wrong_method_warns(run_in_tmpdir, method):
    """Supplying weight windows to a method other than material_wise issues a
    warning and is otherwise ignored."""
    model = _minimal_model()

    # Pre-create the MGXS library so that the conversion skips generation
    # entirely (no simulation is run, so no nuclear data is required).
    Path("mgxs.h5").touch()

    with pytest.warns(UserWarning, match="weight_windows_file"):
        model.convert_to_multigroup(
            method=method,
            weight_windows_file="weight_windows.h5",
        )

    assert model.settings.energy_mode == "multi-group"


def test_weight_windows_file_missing_raises(run_in_tmpdir):
    """A material-wise conversion with a nonexistent weight windows file raises
    a clear error before any simulation is attempted."""
    model = _minimal_model()

    with pytest.raises(FileNotFoundError, match="[Ww]eight windows"):
        model.convert_to_multigroup(
            method="material_wise",
            weight_windows_file="does_not_exist.h5",
        )
