import pytest
import openmc

from tests.testing_harness import PyAPITestHarness


@pytest.fixture
def model():
    model = openmc.model.Model()

    fuel = openmc.Material()
    fuel.set_density('g/cc', 10.0)
    fuel.add_nuclide('U235', 1.0)

    h1 = openmc.Material()
    h1.set_density('g/cc', 0.1)
    h1.add_nuclide('H1', 0.1)

    inner_sphere = openmc.Sphere(x0=1.0, r=5.0)

    fuel_cell = openmc.Cell(fill=fuel, region=-inner_sphere)
    hydrogen_cell = openmc.Cell(fill=h1, region=+inner_sphere)
    univ = openmc.Universe(cells=[fuel_cell, hydrogen_cell])

    # Create one cell on top of the other. Only one
    # has a rotation
    box = openmc.model.RectangularPrism(15., 15., 'z', boundary_type='vacuum')
    lower_z = openmc.ZPlane(-7.5, boundary_type='vacuum')
    upper_z = openmc.ZPlane(22.5, boundary_type='vacuum')
    middle_z = openmc.ZPlane(7.5)

    lower_cell = openmc.Cell(fill=univ, region=-box & +lower_z & -middle_z)
    lower_cell.rotation = (10, 20, 30)
    upper_cell = openmc.Cell(fill=univ, region=-box & +middle_z & -upper_z)
    upper_cell.translation = (0, 0, 15)

    model.geometry = openmc.Geometry(root=[lower_cell, upper_cell])

    model.settings.particles = 10000
    model.settings.inactive = 5
    model.settings.batches = 10
    source_box = openmc.stats.Box((-4., -4., -4.), (4., 4., 4.))
    model.settings.source = openmc.IndependentSource(space=source_box)

    return model


def test_rotation(model):
    harness = PyAPITestHarness('statepoint.10.h5', model)
    harness.main()
