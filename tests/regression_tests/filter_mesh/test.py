import numpy as np
from math import pi,sqrt

import openmc
import pytest

from tests.testing_harness import HashedPyAPITestHarness


@pytest.fixture
def model():
    model = openmc.model.Model()

    fuel = openmc.Material()
    fuel.set_density('g/cm3', 10.0)
    fuel.add_nuclide('U235', 1.0)
    zr = openmc.Material()
    zr.set_density('g/cm3', 1.0)
    zr.add_nuclide('Zr90', 1.0)
    model.materials.extend([fuel, zr])

    box1 = openmc.model.RectangularPrism(10.0, 10.0)
    box2 = openmc.model.RectangularPrism(20.0, 20.0, boundary_type='reflective')
    top = openmc.ZPlane(z0=10.0, boundary_type='vacuum')
    bottom = openmc.ZPlane(z0=-10.0, boundary_type='vacuum')
    cell1 = openmc.Cell(fill=fuel, region=-box1 & +bottom & -top)
    cell2 = openmc.Cell(fill=zr, region=+box1 & -box2 & +bottom & -top)
    model.geometry = openmc.Geometry([cell1, cell2])

    model.settings.batches = 5
    model.settings.inactive = 0
    model.settings.particles = 1000

    # Create meshes
    mesh_1d = openmc.RegularMesh()
    mesh_1d.dimension = [5]
    mesh_1d.lower_left = [-7.5]
    mesh_1d.upper_right = [7.5]

    mesh_2d = openmc.RegularMesh()
    mesh_2d.dimension = [5, 5]
    mesh_2d.lower_left = [-7.5, -7.5]
    mesh_2d.upper_right = [7.5, 7.5]

    mesh_3d = openmc.RegularMesh()
    mesh_3d.dimension = [5, 5, 5]
    mesh_3d.lower_left = [-7.5, -7.5, -7.5]
    mesh_3d.upper_right = [7.5, 7.5, 7.5]
    dx = dy = dz = 15 / 5
    reg_mesh_exp_vols = np.full(mesh_3d.dimension, dx*dy*dz)
    np.testing.assert_equal(mesh_3d.volumes, reg_mesh_exp_vols)

    recti_mesh = openmc.RectilinearMesh()
    recti_mesh.x_grid = np.linspace(-7.5, 7.5, 18)
    recti_mesh.y_grid = np.linspace(-7.5, 7.5, 18)
    recti_mesh.z_grid = np.logspace(0, np.log10(7.5), 11)
    dx = dy = 15 / 17
    dz = np.diff(np.logspace(0, np.log10(7.5), 11))
    dxdy = np.full(recti_mesh.dimension[:2], dx*dy)
    recti_mesh_exp_vols = np.multiply.outer(dxdy, dz)
    np.testing.assert_allclose(recti_mesh.volumes, recti_mesh_exp_vols)

    cyl_mesh = openmc.CylindricalMesh(
        r_grid=np.linspace(0, 7.5, 18),
        phi_grid=np.linspace(0, 2*pi, 19),
        z_grid=np.linspace(-7.5, 7.5, 17),
    )
    dr = 0.5 * np.diff(np.linspace(0, 7.5, 18)**2)
    dp = np.full(cyl_mesh.dimension[1], 2*pi / 18)
    dz = np.full(cyl_mesh.dimension[2], 15 / 16)
    drdp = np.outer(dr, dp)
    cyl_mesh_exp_vols = np.multiply.outer(drdp, dz)
    np.testing.assert_allclose(cyl_mesh.volumes, cyl_mesh_exp_vols)

    sph_mesh = openmc.SphericalMesh(
        r_grid=np.linspace(0, 7.5, 18),
        theta_grid=np.linspace(0, pi, 9),
        phi_grid=np.linspace(0, 2*pi, 19)
    )
    dr = np.diff(np.linspace(0, 7.5, 18)**3) / 3
    dt = np.diff(-np.cos(np.linspace(0, pi, 9)))
    dp = np.full(sph_mesh.dimension[2], 2*pi / 18)
    drdt = np.outer(dr, dt)
    sph_mesh_exp_vols = np.multiply.outer(drdt, dp)
    np.testing.assert_allclose(sph_mesh.volumes, sph_mesh_exp_vols)

    hex_mesh = openmc.HexagonalMesh()
    hex_mesh.dimension = [5,1]
    hex_mesh.lower_left = [-7.5, -7.5]
    hex_mesh.upper_right = [7.5, 7.5]
    dx = 15 / 5
    dz = 15
    dxdy = dx * dx * sqrt(3)
    hm_mesh_exp_vols = dxdy * dz
    np.testing.assert_allclose(hex_mesh.volumes, hm_mesh_exp_vols)

    # Create filters
    reg_filters = [
        openmc.MeshFilter(mesh_1d),
        openmc.MeshFilter(mesh_2d),
        openmc.MeshFilter(mesh_3d),
        openmc.MeshFilter(recti_mesh),
        openmc.MeshFilter(cyl_mesh),
        openmc.MeshFilter(sph_mesh),
        openmc.MeshFilter(hex_mesh),
    ]
    surf_filters = [
        openmc.MeshSurfaceFilter(mesh_1d),
        openmc.MeshSurfaceFilter(mesh_2d),
        openmc.MeshSurfaceFilter(mesh_3d),
        openmc.MeshSurfaceFilter(recti_mesh),
        openmc.MeshSurfaceFilter(cyl_mesh),
        openmc.MeshSurfaceFilter(sph_mesh),
        openmc.MeshSurfaceFilter(hex_mesh),
    ]

    # Create tallies
    for f1, f2 in zip(reg_filters, surf_filters):
        tally = openmc.Tally()
        tally.filters = [f1]
        tally.scores = ['total']
        model.tallies.append(tally)
        tally = openmc.Tally()
        tally.filters = [f2]
        tally.scores = ['current']
        model.tallies.append(tally)

    return model


def test_filter_mesh(model):
    harness = HashedPyAPITestHarness('statepoint.5.h5', model)
    harness.main()
