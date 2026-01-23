import pytest
import openmc

from tests.testing_harness import PyAPITestHarness


@pytest.fixture
def water_model():
    # Define materals and geometry
    water = openmc.Material()
    water.add_nuclide("H1", 2.0)
    water.add_nuclide("O16", 1.0)
    water.set_density("g/cc", 1.0)
    sphere = openmc.Sphere(r=1.0, boundary_type="reflective")
    sph = openmc.Cell(fill=water, region=-sphere)
    geometry = openmc.Geometry([sph])
    source = openmc.IndependentSource(
        energy=openmc.stats.delta_function(10.0e6),
        particle="electron"
    )

    # Define settings
    settings = openmc.Settings()
    settings.particles = 10000
    settings.batches = 1
    settings.cutoff = {"energy_photon": 1000.0}
    settings.run_mode = "fixed source"
    settings.source = source

    # Define tallies
    tally = openmc.Tally()
    tally.scores = ["heating"]
    tallies = openmc.Tallies([tally])

    return openmc.Model(geometry=geometry, settings=settings, tallies=tallies)


def test_electron_heating_calc(water_model):
    harness = PyAPITestHarness("statepoint.1.h5", water_model)
    harness.main()
