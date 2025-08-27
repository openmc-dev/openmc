import numpy as np
import openmc
from openmc.utility_funcs import change_directory
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
    
    lattice00 = openmc.RectLattice()
    lattice00.lower_left = (-d, -d)
    lattice00.pitch = (d, d)
    lattice00.outer = pin
    lattice00.universes = [[pin]]
    box00 = openmc.model.RectangularPrism(d, d, origin=(-d/2,-d/2))
    
    lattice01 = openmc.RectLattice()
    lattice01.lower_left = (-d, 0)
    lattice01.pitch = (d, d)
    lattice01.outer = pin
    lattice01.universes = [[pin]]
    box01 = openmc.model.RectangularPrism(d, d, origin=(-d/2,d/2))
    
    lattice10 = openmc.RectLattice()
    lattice10.lower_left = (0, -d)
    lattice10.pitch = (d, d)
    lattice10.outer = pin
    lattice10.universes = [[pin]]
    box10 = openmc.model.RectangularPrism(d, d, origin=(d/2,-d/2))
    
    lattice11 = openmc.RectLattice()
    lattice11.lower_left = (0, 0)
    lattice11.pitch = (d, d)
    lattice11.outer = pin
    lattice11.universes = [[pin]]
    box11 = openmc.model.RectangularPrism(d, d, origin=(d/2,d/2))
    
    
    cell00 = openmc.Cell(fill=lattice00, region = -box00)
    cell01 = openmc.Cell(fill=lattice01, region = -box01)
    cell10 = openmc.Cell(fill=lattice10, region = -box10)
    cell11 = openmc.Cell(fill=lattice11, region = -box11)
    
    univ = openmc.Universe(cells=[cell00, cell01, cell10, cell11])
    
    box = openmc.model.RectangularPrism(2*d, 2*d, boundary_type='reflective')

    main_cell = openmc.Cell(fill=univ, region=-box)
    model.geometry = openmc.Geometry([main_cell])
    model.geometry.merge_surfaces = True
   
    model.settings.batches = 10
    model.settings.inactive = 5
    model.settings.particles = 1000
    
    return model

@pytest.mark.parametrize("distribmat", [False, True])
def test_lattice(model, distribmat):
    with change_directory(str(distribmat)):
        openmc.reset_auto_ids()
        if distribmat:
            model.differentiate_mats(depletable_only=False)
        harness = PyAPITestHarness('statepoint.10.h5', model)
        harness.main()
