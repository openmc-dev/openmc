import openmc
import pytest

from tests.testing_harness import PyAPITestHarness


@pytest.fixture
def model():
    model = openmc.model.Model()

    m = openmc.Material()
    m.set_density('g/cm3', 4.5)
    m.add_nuclide('U235', 1.0)
    model.materials.append(m)

    sph = openmc.Sphere(r=10.0, boundary_type='vacuum')
    c = openmc.Cell(fill=m, region=-sph)
    model.geometry = openmc.Geometry([c])

    model.settings.particles = 1000
    model.settings.inactive = 3
    model.settings.batches = 7
    model.settings.generations_per_batch = 3
    space = openmc.stats.Box((-4.0, -4.0, -4.0), (4.0, 4.0, 4.))
    model.settings.source = openmc.Source(space=space)

    t = openmc.Tally()
    t.scores = ['flux']
    model.tallies.append(t)

    return model


def test_eigenvalue_genperbatch(model):
    harness = PyAPITestHarness('statepoint.7.h5', model)
    harness.main()
