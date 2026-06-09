"""Tests for the ``weight_windows_file`` argument of
:meth:`openmc.Model.convert_to_multigroup`.

These tests focus on the input validation and the plumbing that threads a weight
windows file into the material-wise MGXS generation simulation. They are
deliberately lightweight and do not run a transport simulation (and therefore do
not require nuclear data): the cross section generation step is either skipped
(by pre-creating the library file) or intercepted via monkeypatching.
"""

from pathlib import Path

import pytest

import openmc


def _sphere_model():
    """Build a small, self-contained fixed source model."""
    openmc.reset_auto_ids()
    model = openmc.Model()

    mat = openmc.Material(name="fuel")
    mat.add_nuclide("U235", 1.0)
    mat.set_density("g/cm3", 1.0)
    model.materials.append(mat)

    sph = openmc.Sphere(r=10.0, boundary_type="vacuum")
    cell = openmc.Cell(fill=mat, region=-sph)
    model.geometry = openmc.Geometry([cell])

    model.settings.run_mode = "fixed source"
    model.settings.particles = 100
    model.settings.batches = 10
    model.settings.source = openmc.IndependentSource(space=openmc.stats.Point())
    return model


def test_weight_windows_file_threaded_into_material_wise(run_in_tmpdir, monkeypatch):
    """The provided weight windows file is loaded and enabled on the model that
    is actually run during material-wise MGXS generation, with its path resolved
    to an absolute location (the solve runs in a temporary directory)."""
    model = _sphere_model()

    # A file must exist for the material_wise validation to pass.
    Path("weight_windows.h5").touch()

    captured = {}

    class _StopGeneration(Exception):
        pass

    def spy(run_model, groups, correction, directory):
        captured["file"] = run_model.settings.weight_windows_file
        captured["on"] = run_model.settings.weight_windows_on
        # Short-circuit before any transport simulation is launched.
        raise _StopGeneration

    monkeypatch.setattr(
        openmc.Model, "_auto_generate_mgxs_lib", staticmethod(spy))

    with pytest.raises(_StopGeneration):
        model.convert_to_multigroup(
            method="material_wise",
            groups="CASMO-2",
            weight_windows_file="weight_windows.h5",
            overwrite_mgxs_library=True,
        )

    assert captured["on"] is True
    assert captured["file"] is not None
    ww_path = Path(captured["file"])
    assert ww_path.is_absolute()
    assert ww_path.name == "weight_windows.h5"


@pytest.mark.parametrize("method", ["stochastic_slab", "infinite_medium"])
def test_weight_windows_file_wrong_method_warns(run_in_tmpdir, method):
    """Supplying weight windows to a method other than material_wise issues a
    warning and is otherwise ignored."""
    model = _sphere_model()

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
    """A material-wise conversion with a nonexistent weight windows file raises a
    clear error before any simulation is attempted."""
    model = _sphere_model()

    with pytest.raises(FileNotFoundError, match="[Ww]eight windows"):
        model.convert_to_multigroup(
            method="material_wise",
            weight_windows_file="does_not_exist.h5",
        )
