import openmc
from pytest import fixture

from tests.testing_harness import PyAPITestHarness

@fixture
def heating_model():
    mat = openmc.Material()
    mat.add_nuclide("Cu63", 0.5)
    mat.add_nuclide("Cu65", 0.5)
    mat.set_density("g/cm3", 1.0)

    sphere = openmc.Sphere(r=20, boundary_type="reflective")
    inside_sphere = openmc.Cell(fill=mat, region=-sphere)
    model = openmc.Model()
    model.geometry = openmc.Geometry([inside_sphere])

    model.settings.particles = 10000
    model.settings.batches = 1
    model.settings.photon_transport = True
    model.settings.electron_treatment = "ttb"
    model.settings.cutoff = {"energy_photon": 1000}
    model.settings.run_mode = "fixed source"
    model.settings.source = openmc.IndependentSource(
        energy=openmc.stats.delta_function(10.0e6),
        particle="photon"
    )

    tally1 = openmc.Tally()
    tally1.scores = ["heating"]
    tally1.nuclides = ["Cu63", "Cu65"]
    tally2 = openmc.Tally()
    tally2.scores = ["heating"]
    model.tallies = [tally1, tally2]

    return model


def test_nuclide_heating(heating_model):
    harness = PyAPITestHarness("statepoint.1.h5", model=heating_model)
    harness.main()
