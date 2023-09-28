from math import pi

import numpy as np
import pytest
import openmc


@pytest.mark.parametrize("val_left,val_right", [(0, 0), (-1., -1.), (2.0, 2)])
def test_raises_error_when_flat(val_left, val_right):
    """Checks that an error is raised when a mesh is flat"""
    mesh = openmc.RegularMesh()

    # Same X
    with pytest.raises(ValueError):
        mesh.lower_left = [val_left, -25, -25]
        mesh.upper_right = [val_right, 25, 25]

    with pytest.raises(ValueError):
        mesh.upper_right = [val_right, 25, 25]
        mesh.lower_left = [val_left, -25, -25]

    # Same Y
    with pytest.raises(ValueError):
        mesh.lower_left = [-25, val_left, -25]
        mesh.upper_right = [25, val_right, 25]

    with pytest.raises(ValueError):
        mesh.upper_right = [25, val_right, 25]
        mesh.lower_left = [-25, val_left, -25]

    # Same Z
    with pytest.raises(ValueError):
        mesh.lower_left = [-25, -25, val_left]
        mesh.upper_right = [25, 25, val_right]

    with pytest.raises(ValueError):
        mesh.upper_right = [25, 25, val_right]
        mesh.lower_left = [-25, -25, val_left]


def test_regular_mesh_bounding_box():
    mesh = openmc.RegularMesh()
    mesh.lower_left = [-2, -3, -5]
    mesh.upper_right = [2, 3, 5]
    bb = mesh.bounding_box
    assert isinstance(bb, openmc.BoundingBox)
    np.testing.assert_array_equal(bb.lower_left, (-2, -3 ,-5))
    np.testing.assert_array_equal(bb.upper_right, (2, 3, 5))


def test_cylindrical_mesh_bounding_box():
    # test with mesh at origin (0, 0, 0)
    mesh = openmc.CylindricalMesh(
        r_grid=[0.1, 0.2, 0.5, 1.],
        z_grid=[0.1, 0.2, 0.4, 0.6, 1.],
        origin=(0, 0, 0)
    )
    np.testing.assert_array_equal(mesh.upper_right, (1, 1, 1))
    np.testing.assert_array_equal(mesh.lower_left, (-1, -1, 0.1))
    bb = mesh.bounding_box
    assert isinstance(bb, openmc.BoundingBox)
    np.testing.assert_array_equal(bb.lower_left, (-1, -1, 0.1))
    np.testing.assert_array_equal(bb.upper_right, (1, 1, 1))

    # test with mesh at origin (3, 5, 7)
    mesh.origin = (3, 5, 7)
    np.testing.assert_array_equal(mesh.upper_right, (4, 6, 8))
    np.testing.assert_array_equal(mesh.lower_left, (2, 4, 7.1))
    bb = mesh.bounding_box
    assert isinstance(bb, openmc.BoundingBox)
    np.testing.assert_array_equal(bb.lower_left, (2, 4, 7.1))
    np.testing.assert_array_equal(bb.upper_right, (4, 6, 8))

    # changing z grid to contain negative numbers
    mesh.z_grid = [-10, 0, 10]
    np.testing.assert_array_equal(mesh.lower_left, (2, 4, -3))
    np.testing.assert_array_equal(mesh.upper_right, (4, 6, 17))


def test_spherical_mesh_bounding_box():
    # test with mesh at origin (0, 0, 0)
    mesh = openmc.SphericalMesh([0.1, 0.2, 0.5, 1.], origin=(0., 0., 0.))
    np.testing.assert_array_equal(mesh.upper_right, (1, 1, 1))
    np.testing.assert_array_equal(mesh.lower_left, (-1, -1, -1))
    bb = mesh.bounding_box
    assert isinstance(bb, openmc.BoundingBox)
    np.testing.assert_array_equal(bb.lower_left, (-1, -1, -1))
    np.testing.assert_array_equal(bb.upper_right, (1, 1, 1))

    # test with mesh at origin (3, 5, 7)
    mesh.origin = (3, 5, 7)
    np.testing.assert_array_equal(mesh.upper_right, (4, 6, 8))
    np.testing.assert_array_equal(mesh.lower_left, (2, 4, 6))
    bb = mesh.bounding_box
    assert isinstance(bb, openmc.BoundingBox)
    np.testing.assert_array_equal(bb.lower_left, (2, 4, 6))
    np.testing.assert_array_equal(bb.upper_right, (4, 6, 8))


def test_SphericalMesh_initiation():
    # test defaults
    mesh = openmc.SphericalMesh(r_grid=(0, 10))
    assert (mesh.origin == np.array([0, 0, 0])).all()
    assert (mesh.r_grid == np.array([0, 10])).all()
    assert (mesh.theta_grid == np.array([0, pi])).all()
    assert (mesh.phi_grid == np.array([0, 2*pi])).all()

    # test setting on creation
    mesh = openmc.SphericalMesh(
        origin=(1, 2, 3),
        r_grid=(0, 2),
        theta_grid=(1, 3),
        phi_grid=(2, 4)
    )
    assert (mesh.origin == np.array([1, 2, 3])).all()
    assert (mesh.r_grid == np.array([0., 2.])).all()
    assert (mesh.theta_grid == np.array([1, 3])).all()
    assert (mesh.phi_grid == np.array([2, 4])).all()

    # test attribute changing
    mesh.r_grid = (0, 11)
    assert (mesh.r_grid == np.array([0., 11.])).all()

    # waffles and pancakes are unfortunately not valid radii
    with pytest.raises(TypeError):
        openmc.SphericalMesh(('🧇', '🥞'))


def test_CylindricalMesh_initiation():
    # test defaults
    mesh = openmc.CylindricalMesh(r_grid=(0, 10), z_grid=(0, 10))
    assert (mesh.origin == np.array([0, 0, 0])).all()
    assert (mesh.r_grid == np.array([0, 10])).all()
    assert (mesh.phi_grid == np.array([0, 2*pi])).all()
    assert (mesh.z_grid == np.array([0, 10])).all()

    # test setting on creation
    mesh = openmc.CylindricalMesh(
        origin=(1, 2, 3),
        r_grid=(0, 2),
        z_grid=(1, 3),
        phi_grid=(2, 4)
    )
    assert (mesh.origin == np.array([1, 2, 3])).all()
    assert (mesh.r_grid == np.array([0., 2.])).all()
    assert (mesh.z_grid == np.array([1, 3])).all()
    assert (mesh.phi_grid == np.array([2, 4])).all()

    # test attribute changing
    mesh.r_grid = (0., 10.)
    assert (mesh.r_grid == np.array([0, 10.])).all()
    mesh.z_grid = (0., 4.)
    assert (mesh.z_grid == np.array([0, 4.])).all()

    # waffles and pancakes are unfortunately not valid radii
    with pytest.raises(TypeError):
        openmc.SphericalMesh(('🧇', '🥞'))
