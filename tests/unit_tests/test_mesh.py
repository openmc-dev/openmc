from math import pi

import numpy as np
import pytest
import openmc

from openmc import RegularMesh, RectilinearMesh, CylindricalMesh, SphericalMesh

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


def test_centroids():
    # regular mesh
    mesh = openmc.RegularMesh()
    mesh.lower_left = (1., 2., 3.)
    mesh.upper_right = (11., 12., 13.)
    mesh.dimension = (1, 1, 1)
    np.testing.assert_array_almost_equal(mesh.centroids[0, 0, 0], [6., 7., 8.])

    # rectilinear mesh
    mesh = openmc.RectilinearMesh()
    mesh.x_grid = [1., 11.]
    mesh.y_grid = [2., 12.]
    mesh.z_grid = [3., 13.]
    np.testing.assert_array_almost_equal(mesh.centroids[0, 0, 0], [6., 7., 8.])

    # cylindrical mesh
    mesh = openmc.CylindricalMesh(r_grid=(0, 10), z_grid=(0, 10), phi_grid=(0, np.pi))
    np.testing.assert_array_almost_equal(mesh.centroids[0, 0, 0], [0.0, 5.0, 5.0])
    # ensure that setting an origin is handled correctly
    mesh.origin = (5.0, 0, -10)
    np.testing.assert_array_almost_equal(mesh.centroids[0, 0, 0], [5.0, 5.0, -5.0])

    # spherical mesh, single element xyz-positive octant
    mesh = openmc.SphericalMesh(r_grid=[0, 10], theta_grid=[0, 0.5*np.pi], phi_grid=[0, np.pi])
    x = 5.*np.cos(0.5*np.pi)*np.sin(0.25*np.pi)
    y = 5.*np.sin(0.5*np.pi)*np.sin(0.25*np.pi)
    z = 5.*np.sin(0.25*np.pi)
    np.testing.assert_array_almost_equal(mesh.centroids[0, 0, 0], [x, y, z])

    mesh.origin = (-5.0, -5.0, 5.0)
    np.testing.assert_array_almost_equal(mesh.centroids[0, 0, 0], [x-5.0, y-5.0, z+5.0])


@pytest.mark.parametrize('mesh_type', ('regular', 'rectilinear', 'cylindrical', 'spherical'))
def test_mesh_vertices(mesh_type):

    ijk = (2, 3, 2)

    # create a new mesh object
    if mesh_type == 'regular':
        mesh = openmc.RegularMesh()
        ll = np.asarray([0.]*3)
        width = np.asarray([0.5]*3)
        mesh.lower_left = ll
        mesh.width = width
        mesh.dimension = (5, 7, 9)

        # spot check that an element has the correct vertex coordinates asociated with it
        # (using zero-indexing here)
        exp_i_j_k = ll + np.asarray(ijk, dtype=float) * width
        np.testing.assert_equal(mesh.vertices[ijk], exp_i_j_k)

        # shift the mesh using the llc
        shift  = np.asarray((3.0, 6.0, 10.0))
        mesh.lower_left += shift
        np.testing.assert_equal(mesh.vertices[ijk], exp_i_j_k+shift)
    elif mesh_type == 'rectilinear':
        mesh = openmc.RectilinearMesh()
        w = np.asarray([0.5] * 3)
        ll = np.asarray([0.]*3)
        dims = (5, 7, 9)
        mesh.x_grid = np.linspace(ll[0], w[0]*dims[0], dims[0])
        mesh.y_grid = np.linspace(ll[1], w[1]*dims[1], dims[1])
        mesh.z_grid = np.linspace(ll[2], w[2]*dims[2], dims[2])
        exp_vert = np.asarray((mesh.x_grid[2], mesh.y_grid[3], mesh.z_grid[2]))
        np.testing.assert_equal(mesh.vertices[ijk], exp_vert)
    elif mesh_type == 'cylindrical':
        r_grid = np.linspace(0, 5, 10)
        z_grid = np.linspace(-10, 10, 20)
        phi_grid = np.linspace(0, 2*np.pi, 8)
        mesh = openmc.CylindricalMesh(r_grid=r_grid, z_grid=z_grid, phi_grid=phi_grid)
        exp_vert = np.asarray((mesh.r_grid[2], mesh.phi_grid[3], mesh.z_grid[2]))
        np.testing.assert_equal(mesh.vertices_cylindrical[ijk], exp_vert)
    elif mesh_type == 'spherical':
        r_grid = np.linspace(0, 13, 14)
        theta_grid = np.linspace(0, np.pi, 11)
        phi_grid = np.linspace(0, 2*np.pi, 7)
        mesh = openmc.SphericalMesh(r_grid=r_grid, theta_grid=theta_grid, phi_grid=phi_grid)
        exp_vert = np.asarray((mesh.r_grid[2], mesh.theta_grid[3], mesh.phi_grid[2]))
        np.testing.assert_equal(mesh.vertices_spherical[ijk], exp_vert)



