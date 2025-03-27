import pytest
import openmc

from tests.testing_harness import PyAPITestHarness


@pytest.fixture
def model():
    u = openmc.Material()
    u.add_nuclide('U235', 1.0)
    u.set_density('g/cm3', 4.5)
    sph = openmc.Sphere(r=10.0, boundary_type='vacuum')
    cell = openmc.Cell(fill=u, region=-sph)
    model = openmc.Model()
    model.geometry = openmc.Geometry([cell])

    model.settings.batches = 10
    model.settings.inactive = 5
    model.settings.particles = 1000
    model.settings.stride = 1_529_170
    model.settings.source = openmc.IndependentSource(
        space=openmc.stats.Box([-4, -4, -4], [4, 4, 4])
    )
    return model


def test_seed(model):
    harness = PyAPITestHarness('statepoint.10.h5', model)
    harness.main()
