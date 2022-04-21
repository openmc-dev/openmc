import numpy as np
import openmc
import pytest

from tests.testing_harness import PyAPITestHarness


def uniform_stack_lattice_model():
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

    rc = 0.4
    h = 1.5
    cyl = openmc.ZCylinder(r=rc)
    top = openmc.ZPlane(z0=h)
    bottom = openmc.ZPlane(z0=0.)
    pellet = -cyl & -top & +bottom
    water_slice = +cyl & -top & +bottom

    fuel = openmc.Cell(fill=uo2, region=pellet)
    water_reflector = openmc.Cell(fill=water, region=water_slice)
    layer = openmc.Universe(cells=[fuel, water_reflector])

    n_pellets = 200

    top = openmc.ZPlane(z0=n_pellets * h)
    tb_refl = openmc.Cell(fill=water, region=-bottom | +top)

    univs = [layer] * n_pellets
    pellet_stack = openmc.StackLattice()
    pellet_stack.central_axis = (0., 0.)
    pellet_stack.base_coordinate = 0.
    pellet_stack.universes = univs
    pellet_stack.is_uniform = True
    pellet_stack.pitch = h

    stack_cell = openmc.Cell(fill=pellet_stack)

    pin_univ = openmc.Universe(cells=[stack_cell, tb_refl])

    d = 1.5 * rc
    box = openmc.model.RectangularParallelepiped(-d, d, -d, d,
                                                 0. - d, n_pellets * h + d,
                                                 boundary_type='reflective')

    main_cell = openmc.Cell(fill=pin_univ, region=-box)
    model.geometry = openmc.Geometry([main_cell])

    model.settings.batches = 10
    model.settings.inactive = 5
    model.settings.particles = 1000

    return model


def test_lattice_stack_uniform():
    model = uniform_stack_lattice_model()
    harness = PyAPITestHarness('statepoint.10.h5', model)
    harness.main()
