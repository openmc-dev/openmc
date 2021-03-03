import numpy as np
import openmc

from tests.testing_harness import PyAPITestHarness


def rotated_lattice_model():
    model = openmc.model.Model()

    # Create some materials
    fuel1 = openmc.Material()
    fuel1.set_density('g/cm3', 10.0)
    fuel1.add_nuclide('U235', 1.0)
    fuel2 = openmc.Material()
    fuel2.set_density('g/cm3', 10.0)
    fuel2.add_nuclide('U238', 1.0)
    water = openmc.Material()
    water.set_density('g/cm3', 1.0)
    water.add_nuclide('H1', 2.0)
    water.add_nuclide('O16', 1.0)
    water.add_s_alpha_beta('c_H_in_H2O')
    model.materials.extend([fuel1, fuel2, water])

    # Create universes for lattices
    r_pin = openmc.ZCylinder(r=0.25)
    fuel_cell = openmc.Cell(fill=fuel1, region=-r_pin)
    water_cell = openmc.Cell(fill=water, region=+r_pin)
    pin_universe = openmc.Universe(cells=(fuel_cell, water_cell))
    r_big_pin = openmc.ZCylinder(r=0.5)
    fuel2_cell = openmc.Cell(fill=fuel2, region=-r_big_pin)
    water2_cell = openmc.Cell(fill=water, region=+r_big_pin)
    big_pin_universe = openmc.Universe(cells=(fuel2_cell, water2_cell))
    all_water_cell = openmc.Cell(fill=water)
    outer_universe = openmc.Universe(30, cells=(all_water_cell,))

    # Create hexagonal lattice
    pitch = 1.25
    hexlat = openmc.HexLattice()
    hexlat.center = (0., 0.)
    hexlat.pitch = [pitch]
    hexlat.outer = outer_universe
    outer_ring = [big_pin_universe] + [pin_universe]*11
    middle_ring = [big_pin_universe] + [pin_universe]*5
    inner_ring = [big_pin_universe]
    hexlat.universes = [outer_ring, middle_ring, inner_ring]

    # Create rectangular lattice
    rectlat = openmc.RectLattice()
    rectlat.center = (0., 0.)
    rectlat.pitch = (pitch, pitch)
    rectlat.lower_left = (-2*pitch, -2*pitch)
    rectlat.outer = outer_universe
    rectlat.universes = np.full((4, 4), pin_universe)
    rectlat.universes[0] = big_pin_universe

    # Create cell filled with translated/rotated rectangular lattice on left
    left_cyl = openmc.ZCylinder(x0=-4.0, r=4.0)
    left_cell = openmc.Cell(fill=rectlat, region=-left_cyl)
    left_cell.translation = (-4.0, 0.0, 0.0)
    left_cell.rotation = (0.0, 0.0, 45.0)

    # Create cell filled with translated/rotated hexagonal lattice on right
    right_cyl = openmc.ZCylinder(x0=4.0, r=4.0)
    right_cell = openmc.Cell(fill=hexlat, region=-right_cyl)
    right_cell.translation = (4.0, 0.0, 0.0)
    right_cell.rotation = (0.0, 0.0, 30.0)

    # Finish up with the geometry
    outer_cyl = openmc.ZCylinder(r=8.0, boundary_type='vacuum')
    main_cell = openmc.Cell(fill=water, region=-outer_cyl & +left_cyl & +right_cyl)
    model.geometry = openmc.Geometry([main_cell, left_cell, right_cell])

    model.settings.batches = 5
    model.settings.inactive = 0
    model.settings.particles = 1000
    model.settings.source = openmc.Source(space=openmc.stats.Point())
    model.settings.export_to_xml()

    return model


def test():
    model = rotated_lattice_model()
    harness = PyAPITestHarness('statepoint.5.h5', model)
    harness.main()
