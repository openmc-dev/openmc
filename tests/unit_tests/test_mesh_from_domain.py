import numpy as np
import openmc
import pytest


def test_mesh_from_cell():
    """Tests a RegularMesh can be made from a Cell and the specified dimensions
    are propagated through. Cell is not centralized"""
    surface = openmc.Sphere(r=10, x0=2, y0=3, z0=5)
    cell = openmc.Cell(region=-surface)

    mesh = openmc.RegularMesh.from_domain(cell, dimension=[7, 11, 13])
    assert isinstance(mesh, openmc.RegularMesh)
    assert np.array_equal(mesh.dimension, (7, 11, 13))
    assert np.array_equal(mesh.lower_left, cell.bounding_box[0])
    assert np.array_equal(mesh.upper_right, cell.bounding_box[1])


def test_mesh_from_region():
    """Tests a RegularMesh can be made from a Region and the default dimensions
    are propagated through. Region is not centralized"""
    surface = openmc.Sphere(r=1, x0=-5, y0=-3, z0=-2)
    region = -surface

    mesh = openmc.RegularMesh.from_domain(region)
    assert isinstance(mesh, openmc.RegularMesh)
    assert np.array_equal(mesh.dimension, (100, 100, 100))  # default values
    assert np.array_equal(mesh.lower_left, region.bounding_box[0])
    assert np.array_equal(mesh.upper_right, region.bounding_box[1])


def test_mesh_from_universe():
    """Tests a RegularMesh can be made from a Universe and the default dimensions
    are propagated through. Universe is centralized"""
    surface = openmc.Sphere(r=42)
    cell = openmc.Cell(region=-surface)
    universe = openmc.Universe(cells=[cell])

    mesh = openmc.RegularMesh.from_domain(universe)
    assert isinstance(mesh, openmc.RegularMesh)
    assert np.array_equal(mesh.dimension, (100, 100, 100))  # default values
    assert np.array_equal(mesh.lower_left, universe.bounding_box[0])
    assert np.array_equal(mesh.upper_right, universe.bounding_box[1])


def test_mesh_from_geometry():
    """Tests a RegularMesh can be made from a Geometry and the default dimensions
    are propagated through. Geometry is centralized"""
    surface = openmc.Sphere(r=42)
    cell = openmc.Cell(region=-surface)
    universe = openmc.Universe(cells=[cell])
    geometry = openmc.Geometry(universe)

    mesh = openmc.RegularMesh.from_domain(geometry)
    assert isinstance(mesh, openmc.RegularMesh)
    assert np.array_equal(mesh.dimension, (100, 100, 100))  # default values
    assert np.array_equal(mesh.lower_left, geometry.bounding_box[0])
    assert np.array_equal(mesh.upper_right, geometry.bounding_box[1])


def test_error_from_unsupported_object():
    with pytest.raises(TypeError):
        openmc.RegularMesh.from_domain("vacuum energy")
