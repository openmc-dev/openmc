from math import sin, cos, pi

import openmc
import pytest

from tests.testing_harness import PyAPITestHarness


@pytest.fixture
def model():
    model = openmc.model.Model()

    # Define materials
    water = openmc.Material()
    water.add_nuclide('H1', 2.0)
    water.add_nuclide('O16', 1.0)
    water.add_s_alpha_beta('c_H_in_H2O')
    water.set_density('g/cc', 1.0)

    fuel = openmc.Material()
    fuel.add_nuclide('U235', 1.0)
    fuel.set_density('g/cc', 4.5)

    # Define the geometry.  Note that this geometry is somewhat non-sensical
    # (it essentially defines a circle of half-cylinders), but it is
    # designed so that periodic and reflective BCs will give different
    # answers.
    theta1 = (-1/6 + 1/2) * pi
    theta2 = (1/6 - 1/2) * pi
    plane1 = openmc.Plane(a=cos(theta1), b=sin(theta1), boundary_type='periodic')
    plane2 = openmc.Plane(a=cos(theta2), b=sin(theta2), boundary_type='periodic')

    x_max = openmc.XPlane(5., boundary_type='reflective')

    z_cyl = openmc.ZCylinder(x0=3*cos(pi/6), y0=3*sin(pi/6), r=2.0)

    outside_cyl = openmc.Cell(1, fill=water, region=(
        +plane1 & +plane2 & -x_max & +z_cyl))
    inside_cyl = openmc.Cell(2, fill=fuel, region=(
        +plane1 & +plane2 & -z_cyl))
    root_universe = openmc.Universe(0, cells=(outside_cyl, inside_cyl))
    model.geometry = openmc.Geometry(root_universe)

    # Define settings
    model.settings = openmc.Settings()
    model.settings.particles = 1000
    model.settings.batches = 4
    model.settings.inactive = 0
    model.settings.source = openmc.IndependentSource(space=openmc.stats.Box(
        (0, 0, 0), (5, 5, 0))
    )
    return model


def test_periodic(model):
    harness = PyAPITestHarness('statepoint.4.h5', model)
    harness.main()
