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


@pytest.mark.parametrize("num_groups", [None, 6])
def test_get_kinetics_parameters(run_in_tmpdir, num_groups):
    # Create basic model
    model = openmc.Model()
    material = openmc.Material(name="core")
    material.add_nuclide("U235", 1.0)
    material.set_density('g/cm3', 16.0)
    sphere = openmc.Sphere(r=10.0, boundary_type="vacuum")
    cell = openmc.Cell(region=-sphere, fill=material)
    model.geometry = openmc.Geometry([cell])
    model.settings.particles = 1000
    model.settings.batches = 20
    model.settings.inactive = 5
    model.settings.ifp_n_generation = 5

    # Add IFP tallies
    model.add_kinetics_parameters_tallies(num_groups=num_groups)

    # Run and get kinetics parameters
    sp_file = model.run()
    with openmc.StatePoint(sp_file) as sp:
        params = sp.get_kinetics_parameters()
    assert isinstance(params, openmc.KineticsParameters)
    assert params.generation_time is not None
    assert params.beta_effective is not None
    if num_groups is not None:
        assert len(params.beta_effective) == num_groups
