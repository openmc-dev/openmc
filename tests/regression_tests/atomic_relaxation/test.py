import openmc
import pytest

from tests.testing_harness import PyAPITestHarness


@pytest.fixture
def model():
    mat = openmc.Material()
    mat.add_nuclide('Pb208', 1.0)
    mat.set_density('g/cm3', 11.35)

    sphere = openmc.Sphere(r=1.0e9, boundary_type='reflective')
    inside_sphere = openmc.Cell(fill=mat, region=-sphere)
    model = openmc.Model()
    model.geometry = openmc.Geometry([inside_sphere])

    # Isotropic point source of 1 MeV photons at the origin
    model.settings.source = openmc.IndependentSource(
        particle='photon',
        energy=openmc.stats.delta_function(1.0e6)
    )

    # Fixed-source photon transport with atomic relaxation disabled
    model.settings.particles = 10000
    model.settings.batches = 1
    model.settings.photon_transport = True
    model.settings.electron_treatment = 'led'
    model.settings.atomic_relaxation = False
    model.settings.run_mode = 'fixed source'

    tally = openmc.Tally()
    tally.filters = [openmc.ParticleFilter(['photon', 'electron'])]
    tally.scores = ['flux', 'heating']
    model.tallies = [tally]
    return model


def test_atomic_relaxation(model):
    harness = PyAPITestHarness('statepoint.1.h5', model=model)
    harness.main()
