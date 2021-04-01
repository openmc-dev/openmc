import openmc
import pytest

from tests.testing_harness import TestHarness, PyAPITestHarness

@pytest.fixture
def model():
    model = openmc.model.Model()

    fuel = openmc.Material()
    fuel.set_density('g/cc', 10.0)
    fuel.add_nuclide('U235', 1.0)

    h1 = openmc.Material()
    h1.set_density('g/cc', 0.1)
    h1.add_nuclide('H1', 0.1)

    outer_sphere = openmc.Sphere(r=10.0, boundary_type='vacuum')
    inner_sphere = openmc.Sphere(x0=1.0, r=5.0)

    fuel_cell = openmc.Cell(fill=fuel, region=-inner_sphere)
    hydrogen_cell = openmc.Cell(fill=h1, region=+inner_sphere)
    univ = openmc.Universe(cells=[fuel_cell, hydrogen_cell])
    outer_cell = openmc.Cell(fill=univ, region=-outer_sphere)
    outer_cell.rotation = (10, 20, 30)

    model.geometry = openmc.Geometry(root=[outer_cell])

    model.settings.particles = 1000
    model.settings.inactive = 5
    model.settings.batches = 10
    source_box = openmc.stats.Box((-4., -4., -4.), (4., 4., 4.))
    model.settings.source = openmc.Source(space=source_box)

    return model


def test_rotation(model):
    harness = PyAPITestHarness('statepoint.10.h5', model)
    harness.main()
