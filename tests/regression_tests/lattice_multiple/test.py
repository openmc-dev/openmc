import numpy as np
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
    big_cyl = openmc.ZCylinder(r=0.5)
    pin = openmc.model.pin([cyl], [uo2, water])
    big_pin = openmc.model.pin([big_cyl], [uo2, water])

    d = 1.2
    inner_lattice = openmc.RectLattice()
    inner_lattice.lower_left = (-d, -d)
    inner_lattice.pitch = (d, d)
    inner_lattice.outer = pin
    inner_lattice.universes = [
        [big_pin, pin],
        [pin, pin],
    ]
    inner_cell = openmc.Cell(fill=inner_lattice)
    inner_univ = openmc.Universe(cells=[inner_cell])

    lattice = openmc.RectLattice()
    lattice.lower_left = (-2*d, -2*d)
    lattice.pitch = (2*d, 2*d)
    lattice.universes = np.full((2, 2), inner_univ)

    box = openmc.model.rectangular_prism(4*d, 4*d, boundary_type='reflective')
    main_cell = openmc.Cell(fill=lattice, region=box)
    model.geometry = openmc.Geometry([main_cell])

    model.settings.batches = 10
    model.settings.inactive = 5
    model.settings.particles = 1000

    return model


def test_lattice_multiple(model):
    harness = PyAPITestHarness('statepoint.10.h5', model)
    harness.main()
