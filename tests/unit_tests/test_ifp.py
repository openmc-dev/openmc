"""Test the Iterated Fission Probability (IFP) method to compute adjoint-weighted
kinetics parameters using dedicated tallies."""

import pytest
import openmc


def test_xml_serialization(run_in_tmpdir):
    """Check that a simple use case can be written and read in XML."""
    parameter = {"n_generation": 5}
    settings = openmc.Settings()
    settings.iterated_fission_probability = parameter
    settings.export_to_xml()

    read_settings = openmc.Settings.from_xml()
    assert read_settings.iterated_fission_probability == parameter


@pytest.fixture(scope="module")
def geometry():
    openmc.reset_auto_ids()
    material = openmc.Material()
    sphere = openmc.Sphere(r=1.0, boundary_type="vacuum")
    cell = openmc.Cell(region=-sphere, fill=material)
    return openmc.Geometry([cell])


@pytest.mark.parametrize(
    "settings, error",
    [
        ({"iterated_fission_probability": {"n_generation": 0}}, ValueError),
        ({"run_mode": "fixed source", "iterated_fission_probability": {"n_generation": 5}}, RuntimeError),
        ({"n_inactive": 6, "iterated_fission_probability": {"n_generation": 5}}, RuntimeError)
    ],
)
def test_exceptions(settings, error, run_in_tmpdir, geometry):
    """Test settings configuration that should return an error."""
    with pytest.raises(error):
        settings = openmc.Settings(**settings)
        settings.particles = 100
        settings.batches = 10
        model = openmc.Model(geometry=geometry, settings=settings)
        model.run()
