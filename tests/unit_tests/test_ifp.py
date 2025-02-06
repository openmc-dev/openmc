"""Test the Iterated Fission Probability (IFP) method to compute adjoint-weighted
kinetics parameters using dedicated tallies."""

import pytest
import openmc


def test_xml_serialization(run_in_tmpdir):
    """Check that a simple use case can be written and read in XML."""
    parameter = 5
    settings = openmc.Settings()
    settings.ifp_n_generation = parameter
    settings.export_to_xml()

    read_settings = openmc.Settings.from_xml()
    assert read_settings.ifp_n_generation == parameter


@pytest.fixture(scope="module")
def geometry():
    openmc.reset_auto_ids()
    material = openmc.Material()
    material.add_nuclide("U235", 1.0)
    sphere = openmc.Sphere(r=1.0, boundary_type="vacuum")
    cell = openmc.Cell(region=-sphere, fill=material)
    return openmc.Geometry([cell])


@pytest.mark.parametrize(
    "options, error",
    [
        ({"ifp_n_generation": 0}, ValueError),
        ({"ifp_n_generation": -1}, ValueError),
        ({"run_mode": "fixed source"}, RuntimeError),
        ({"inactive": 5, "ifp_n_generation": 6}, RuntimeError),
        ({"inactive": 9}, RuntimeError)
    ],
)
def test_exceptions(options, error, run_in_tmpdir, geometry):
    """Test settings configuration that should return an error."""
    with pytest.raises(error):
        settings = openmc.Settings(**options)
        settings.particles = 100
        settings.batches = 15
        tally = openmc.Tally(name="ifp-scores")
        tally.scores = ["ifp-time-numerator", "ifp-beta-numerator", "ifp-denominator"]
        tallies = openmc.Tallies([tally])
        model = openmc.Model(geometry=geometry, settings=settings, tallies=tallies)
        model.run()
