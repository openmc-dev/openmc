import openmc
import pytest

from tests.testing_harness import PyAPITestHarness


@pytest.fixture
def box_model():
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

    # Define geometry
    x_min = openmc.XPlane(surface_id=1, x0=0., boundary_type='periodic')
    x_max = openmc.XPlane(surface_id=2, x0=5., boundary_type='reflective')

    y_min = openmc.YPlane(surface_id=3, y0=0., boundary_type='periodic')
    y_max = openmc.YPlane(surface_id=4, y0=5., boundary_type='reflective')
    y_min.periodic_surface = x_min

    z_min = openmc.ZPlane(surface_id=5, z0=-5., boundary_type='periodic')
    z_max = openmc.Plane(surface_id=6, a=0, b=0, c=1, d=5.,
                         boundary_type='periodic')
    z_cyl = openmc.ZCylinder(surface_id=7, x0=2.5, y0=0., r=2.0)

    outside_cyl = openmc.Cell(1, fill=water, region=(
        +x_min & -x_max & +y_min & -y_max & +z_min & -z_max & +z_cyl))
    inside_cyl = openmc.Cell(2, fill=fuel, region=(
        +y_min & +z_min & -z_max & -z_cyl))
    root_universe = openmc.Universe(0, cells=(outside_cyl, inside_cyl))
    model.geometry = openmc.Geometry(root_universe)

    # Define settings
    model.settings.particles = 1000
    model.settings.batches = 4
    model.settings.inactive = 0
    model.settings.source = openmc.IndependentSource(space=openmc.stats.Box(
        (0, 0, 0), (5, 5, 0))
    )
    return model


def test_periodic(box_model):
    harness = PyAPITestHarness('statepoint.4.h5', box_model)
    harness.main()
