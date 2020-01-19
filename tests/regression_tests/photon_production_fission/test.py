import openmc
import pytest

from tests.testing_harness import PyAPITestHarness


@pytest.fixture
def model():
    model = openmc.model.Model()
    mat = openmc.Material()
    mat.set_density('g/cm3', 10.0)
    mat.add_nuclide('U235', 1.0)
    model.materials.append(mat)

    sph = openmc.Sphere(r=100.0, boundary_type='reflective')
    cell = openmc.Cell(fill=mat, region=-sph)
    model.geometry = openmc.Geometry([cell])

    model.settings.particles = 1000
    model.settings.batches = 5
    model.settings.inactive = 2
    model.settings.photon_transport = True
    model.settings.source = openmc.Source(space=openmc.stats.Point((0, 0, 0)))

    particle_filter = openmc.ParticleFilter(['neutron', 'photon'])
    tally_tracklength = openmc.Tally()
    tally_tracklength.filters = [particle_filter]
    tally_tracklength.scores = ['fission', 'heating-local']
    tally_tracklength.nuclides = ['U235', 'total']
    tally_tracklength.estimator = 'tracklength'
    tally_collision = openmc.Tally()
    tally_collision.filters = [particle_filter]
    tally_collision.scores = ['fission', 'heating', 'heating-local']
    tally_collision.nuclides = ['U235', 'total']
    tally_collision.estimator = 'collision'
    tally_analog = openmc.Tally()
    tally_analog.filters = [particle_filter]
    tally_analog.scores = ['fission', 'heating', 'heating-local']
    tally_analog.nuclides = ['U235', 'total']
    tally_analog.estimator = 'analog'
    model.tallies.extend([tally_tracklength, tally_collision, tally_analog])

    return model


def test_photon_production_fission(model):
    harness = PyAPITestHarness('statepoint.5.h5', model)
    harness.main()
