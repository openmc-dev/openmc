import openmc
import numpy
import pytest

from tests.testing_harness import PyAPITestHarness


@pytest.fixture
def xcyl_model():
    model = openmc.Model()
    # Define materials
    fuel = openmc.Material()
    fuel.add_nuclide('U235', 0.2)
    fuel.add_nuclide('U238', 0.8)
    fuel.set_density('g/cc', 19.1)
    model.materials = openmc.Materials([fuel])

    # Define geometry
    # finite cylinder
    x_min = openmc.XPlane(x0=0.0, boundary_type='reflective')
    x_max = openmc.XPlane(x0=20.0, boundary_type='reflective')
    x_cyl = openmc.XCylinder(r=20.0,boundary_type='vacuum')
    # slice cylinder for periodic BC
    periodic_bounding_yplane = openmc.YPlane(y0=0, boundary_type='periodic')
    periodic_bounding_plane = openmc.Plane(
        a=0.0, b=-np.sqrt(3) / 3, c=1, boundary_type='periodic',
    )
    sixth_cyl_cell = openmc.Cell(1, fill=fuel, region = 
        +x_min &- x_max & -x_cyl & +periodic_bounding_yplane & +periodic_bounding_plane)
    periodic_bounding_yplane.periodic_surface = periodic_bounding_plane
    periodic_bounding_plane.periodic_surface = periodic_bounding_yplane

    model.geometry = openmc.Geometry([sixth_cyl_cell])


    # Define settings
    model.settings.particles = 1000
    model.settings.batches = 4
    model.settings.inactive = 0
    model.settings.source = openmc.IndependentSource(space=openmc.stats.Box(
        (0, 0, 0), (20, 20, 20))
    )

@pytest.fixture
def ycyl_model():
    model = openmc.Model()
    # Define materials
    fuel = openmc.Material()
    fuel.add_nuclide('U235', 0.2)
    fuel.add_nuclide('U238', 0.8)
    fuel.set_density('g/cc', 19.1)
    model.materials = openmc.Materials([fuel])

    # Define geometry
    # finite cylinder
    y_min = openmc.YPlane(y0=0.0, boundary_type='reflective')
    y_max = openmc.YPlane(y0=20.0, boundary_type='reflective')
    y_cyl = openmc.YCylinder(r=20.0,boundary_type='vacuum')
    # slice cylinder for periodic BC
    periodic_bounding_xplane = openmc.XPlane(x0=0, boundary_type='periodic')
    periodic_bounding_plane = openmc.Plane(
        a=-np.sqrt(3) / 3, b=0.0, c=1, boundary_type='periodic',
    )
    sixth_cyl_cell = openmc.Cell(1, fill=fuel, region = 
        +y_min &- y_max & -y_cyl & +periodic_bounding_xplane & +periodic_bounding_plane)
    periodic_bounding_xplane.periodic_surface = periodic_bounding_plane
    periodic_bounding_plane.periodic_surface = periodic_bounding_xplane
    model.geometry = openmc.Geometry([sixth_cyl_cell])


    # Define settings
    model.settings.particles = 1000
    model.settings.batches = 4
    model.settings.inactive = 0
    model.settings.source = openmc.IndependentSource(space=openmc.stats.Box(
        (0, 0, 0), (20, 20, 20))
    )


@pytest.mark.parametrize('cyl_model', ['xcyl_model','ycyl_model'])
def test_periodic(cyl_model):
    with change_directory(cyl_model):
        harness = PyAPITestHarness('statepoint.4.h5', cyl_model)
        harness.main()
    assert 