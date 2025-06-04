import numpy as np
import openmc
import pytest


def test_reg_mesh_from_cell():
    """Tests a RegularMesh can be made from a Cell and the specified dimensions
    are propagated through. Cell is not centralized"""
    surface = openmc.Sphere(r=10, x0=2, y0=3, z0=5)
    cell = openmc.Cell(region=-surface)

    mesh = openmc.RegularMesh.from_domain(domain=cell, dimension=[7, 11, 13])
    assert isinstance(mesh, openmc.RegularMesh)
    assert np.array_equal(mesh.dimension, (7, 11, 13))
    assert np.array_equal(mesh.lower_left, cell.bounding_box[0])
    assert np.array_equal(mesh.upper_right, cell.bounding_box[1])


def test_cylindrical_mesh_from_cell():
    """Tests a CylindricalMesh can be made from a Cell and the specified
    dimensions are propagated through."""
    # Cell is not centralized on Z axis
    cy_surface = openmc.ZCylinder(r=50)
    z_surface_1 = openmc.ZPlane(z0=40)
    z_surface_2 = openmc.ZPlane(z0=10)
    cell = openmc.Cell(region=-cy_surface & -z_surface_1 & +z_surface_2)
    mesh = openmc.CylindricalMesh.from_domain(domain=cell, dimension=[2, 4, 3])

    assert isinstance(mesh, openmc.CylindricalMesh)
    assert np.array_equal(mesh.dimension, (2, 4, 3))
    assert np.array_equal(mesh.r_grid, [0., 25., 50.])
    assert np.array_equal(mesh.phi_grid, [0., 0.5*np.pi, np.pi, 1.5*np.pi, 2.*np.pi])
    assert np.array_equal(mesh.z_grid, [0., 10., 20., 30.])
    assert np.array_equal(mesh.origin, [0., 0., 10.])

    # Cell is not centralized on Z or X axis
    cy_surface = openmc.ZCylinder(r=50, x0=100)
    cell = openmc.Cell(region=-cy_surface & -z_surface_1 & +z_surface_2)
    mesh = openmc.CylindricalMesh.from_domain(domain=cell, dimension=[1, 1, 1])

    assert isinstance(mesh, openmc.CylindricalMesh)
    assert np.array_equal(mesh.dimension, (1, 1, 1))
    assert np.array_equal(mesh.r_grid, [0., 150.])
    assert np.array_equal(mesh.origin, [100., 0., 10.])

    # Cell is not centralized on Z, X or Y axis
    cy_surface = openmc.ZCylinder(r=50, x0=100, y0=170)
    cell = openmc.Cell(region=-cy_surface & -z_surface_1 & +z_surface_2)
    mesh = openmc.CylindricalMesh.from_domain(domain=cell, dimension=[1, 1, 1])

    assert isinstance(mesh, openmc.CylindricalMesh)
    assert np.array_equal(mesh.r_grid, [0., 220.])
    assert np.array_equal(mesh.origin, [100., 170., 10.])


def test_reg_mesh_from_region():
    """Tests a RegularMesh can be made from a Region and the default dimensions
    are propagated through. Region is not centralized"""
    surface = openmc.Sphere(r=1, x0=-5, y0=-3, z0=-2)
    region = -surface

    mesh = openmc.RegularMesh.from_domain(domain=region)
    assert isinstance(mesh, openmc.RegularMesh)
    assert np.array_equal(mesh.dimension, (10, 10, 10))  # default values
    assert np.array_equal(mesh.lower_left, region.bounding_box[0])
    assert np.array_equal(mesh.upper_right, region.bounding_box[1])


def test_cylindrical_mesh_from_region():
    """Tests a CylindricalMesh can be made from a Region and the specified
    dimensions and phi_grid_bounds are propagated through. Cell is centralized"""
    cy_surface = openmc.ZCylinder(r=6)
    z_surface_1 = openmc.ZPlane(z0=30)
    z_surface_2 = openmc.ZPlane(z0=-30)
    cell = openmc.Cell(region=-cy_surface & -z_surface_1 & +z_surface_2)
    mesh = openmc.CylindricalMesh.from_domain(
        domain=cell,
        dimension=(6, 2, 3),
        phi_grid_bounds=(0., np.pi)
    )

    assert isinstance(mesh, openmc.CylindricalMesh)
    assert np.array_equal(mesh.dimension, (6, 2, 3))
    assert np.array_equal(mesh.r_grid, [0., 1., 2., 3., 4., 5., 6.])
    assert np.array_equal(mesh.phi_grid, [0., 0.5*np.pi, np.pi])
    assert np.array_equal(mesh.z_grid, [0.0, 20., 40., 60])
    assert np.array_equal(mesh.origin, (0.0, 0.0, -30.))


def test_spherical_mesh_from_region():
    """Tests a SphericalMesh can be made from a Region and the specified
    dimensions are propagated through. Cell is not centralized"""
    sphere = openmc.Sphere(r=5, x0=2, y0=3, z0=4)
    region = -sphere

    geometry = openmc.Geometry(openmc.Universe(cells=[openmc.Cell(region=region)]))

    region_mesh = openmc.SphericalMesh.from_domain(
        domain=region, dimension=(4, 3, 4))
    universe_mesh = openmc.SphericalMesh.from_domain(
        domain=geometry.root_universe, dimension=(4, 3, 4))
    geometry_mesh = openmc.SphericalMesh.from_domain(
        domain=geometry, dimension=(4, 3, 4))


    for mesh in (region_mesh, universe_mesh, geometry_mesh):
        assert isinstance(mesh, openmc.SphericalMesh)
        assert np.array_equal(mesh.dimension, (4, 3, 4))
        assert np.array_equal(mesh.r_grid, [0., 1.25, 2.5, 3.75, 5.0])
        assert np.array_equal(mesh.theta_grid, [0., np.pi/3., 2*np.pi/3., np.pi])
        assert np.array_equal(mesh.phi_grid, [0., np.pi/2., np.pi, 3*np.pi/2., 2*np.pi])
        assert np.array_equal(mesh.origin, (2.0, 3.0, 4.0))

        for p in mesh.centroids.reshape(-1, 3):
            assert p in mesh.bounding_box


def test_reg_mesh_from_universe():
    """Tests a RegularMesh can be made from a Universe and the default
    dimensions are propagated through. Universe is centralized"""
    surface = openmc.Sphere(r=42)
    cell = openmc.Cell(region=-surface)
    universe = openmc.Universe(cells=[cell])

    mesh = openmc.RegularMesh.from_domain(universe)
    assert isinstance(mesh, openmc.RegularMesh)
    assert np.array_equal(mesh.dimension, (10, 10, 10))  # default values
    assert np.array_equal(mesh.lower_left, universe.bounding_box[0])
    assert np.array_equal(mesh.upper_right, universe.bounding_box[1])


def test_reg_mesh_from_geometry():
    """Tests a RegularMesh can be made from a Geometry and the default
    dimensions are propagated through. Geometry is centralized"""
    surface = openmc.Sphere(r=42)
    cell = openmc.Cell(region=-surface)
    universe = openmc.Universe(cells=[cell])
    geometry = openmc.Geometry(universe)

    mesh = openmc.RegularMesh.from_domain(geometry)
    assert isinstance(mesh, openmc.RegularMesh)
    assert np.array_equal(mesh.dimension, (10, 10, 10))  # default values
    assert np.array_equal(mesh.lower_left, geometry.bounding_box[0])
    assert np.array_equal(mesh.upper_right, geometry.bounding_box[1])


def test_error_from_unsupported_object():
    with pytest.raises(TypeError):
        openmc.RegularMesh.from_domain("vacuum energy")
