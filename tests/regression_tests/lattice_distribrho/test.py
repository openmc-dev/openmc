import openmc
import pytest

from tests.testing_harness import PyAPITestHarness

@pytest.fixture
def model():
    model = openmc.model.Model()

    uo2 = openmc.Material(name='UO2')
    uo2.set_density('g/cm3', 10.0)
    uo2.add_nuclide('U235', 1.0)
    uo2.add_nuclide('O16', 2.0)
    water = openmc.Material(name='light water')
    water.add_nuclide('H1', 2.0)
    water.add_nuclide('O16', 1.0)
    water.set_density('g/cm3', 1.0)
    water.add_s_alpha_beta('c_H_in_H2O')
    model.materials.extend([uo2, water])

    cyl = openmc.ZCylinder(r=0.4)
    pin = openmc.model.pin([cyl], [uo2, water])
    d = 1.0

    lattice = openmc.RectLattice()
    lattice.lower_left = (-d, -d)
    lattice.pitch = (d, d)
    lattice.universes = [[pin, pin],
                         [pin, pin]]
    box = openmc.model.RectangularPrism(2.0 * d, 2.0 * d, origin=(0.0,0.0), boundary_type='reflective')

    pin.cells[1].density = [10.0, 20.0, 10.0, 20.0]

    model.geometry = openmc.Geometry([openmc.Cell(fill=lattice, region = -box)])
    model.geometry.merge_surfaces = True

    model.settings.batches = 10
    model.settings.inactive = 5
    model.settings.particles = 1000

    return model

def test_lattice_checkerboard(model):
    harness = PyAPITestHarness('statepoint.10.h5', model)
    harness.main()
