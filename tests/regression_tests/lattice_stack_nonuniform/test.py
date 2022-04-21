from tests.testing_harness import PyAPITestHarness

import numpy as np
import openmc
import pytest


def nonuniform_stack_lattice_model():
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
    h1 = 1.0
    h2 = 2.0
    h_avg = (h1 + h2) / 2
    cyl = openmc.ZCylinder(r=rc)
    top1 = openmc.ZPlane(z0=h1)
    top2 = openmc.ZPlane(z0=h2)
    bottom = openmc.ZPlane(z0=0.)
    pellet1 = -cyl & -top1 & +bottom
    pellet2 = -cyl & -top2 & +bottom
    water_slice1 = +cyl & -top1 & +bottom
    water_slice2 = +cyl & -top2 & +bottom

    fuel1 = openmc.Cell(fill=uo2, region=pellet1)
    fuel2 = openmc.Cell(fill=uo2, region=pellet2)
    water_reflector1 = openmc.Cell(fill=water, region=water_slice1)
    water_reflector2 = openmc.Cell(fill=water, region=water_slice2)
    layer1 = openmc.Universe(cells=[fuel1, water_reflector1])
    layer2 = openmc.Universe(cells=[fuel2, water_reflector2])

    n_pellets = 200

    top = openmc.ZPlane(z0=n_pellets * h_avg)
    tb_refl = openmc.Cell(fill=water, region=-bottom | +top)

    d = 1.5 * rc
    box = openmc.model.RectangularParallelepiped(-d, d, -d, d, 0. - d,
                                                 n_pellets * h_avg + d,
                                                 boundary_type='reflective')

    univs = [layer1, layer2] * int(n_pellets / 2)
    pellet_stack = openmc.StackLattice()
    pellet_stack.central_axis = (0., 0.)
    pellet_stack.base_coordinate = 0.
    pellet_stack.universes = univs
    pellet_stack.is_uniform = False
    pellet_stack.pitch = [h1, h2] * int(n_pellets / 2)

    stack_cell = openmc.Cell(fill=pellet_stack)

    pin_univ = openmc.Universe(cells=[stack_cell, tb_refl])

    main_cell = openmc.Cell(fill=pin_univ, region=-box)
    model.geometry = openmc.Geometry([main_cell])

    model.settings.batches = 10
    model.settings.inactive = 5
    model.settings.particles = 1000

    return model


def test_lattice_stack_nonuniform():
    model = nonuniform_stack_lattice_model()
    harness = PyAPITestHarness('statepoint.10.h5', model)
    harness.main()
