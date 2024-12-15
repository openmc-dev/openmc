import openmc
import pytest

from tests.testing_harness import PyAPITestHarness


@pytest.fixture
def time_model():
    model = openmc.Model()
    time_cutoff = 1e-7

    # A single sphere
    s1 = openmc.Sphere(r=200, boundary_type="vacuum")
    sphere = openmc.Cell()
    sphere.region = -s1
    model.geometry = openmc.Geometry([sphere])

    # Set the running parameters
    settings_file = openmc.Settings()
    settings_file.run_mode = "fixed source"
    settings_file.batches = 10
    settings_file.particles = 100
    settings_file.cutoff = {"time_neutron": time_cutoff}
    settings_file.source = openmc.IndependentSource(
        space=openmc.stats.Point(), energy=openmc.stats.Discrete([1e4], [1])
    )
    model.settings = settings_file

    # Tally flux under time cutoff
    tallies = openmc.Tallies()
    tally = openmc.Tally()
    tally.scores = ["flux"]
    time_filter = openmc.TimeFilter([0, time_cutoff, 2 * time_cutoff])
    tally.filters = [time_filter]
    tallies.append(tally)
    model.tallies = tallies

    return model


def test_time_cutoff(time_model):
    harness = PyAPITestHarness("statepoint.10.h5", time_model)
    harness.main()
