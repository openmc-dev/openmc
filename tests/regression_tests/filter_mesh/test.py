import numpy as np

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

    box1 = openmc.model.rectangular_prism(10.0, 10.0)
    box2 = openmc.model.rectangular_prism(20.0, 20.0, boundary_type='reflective')
    top = openmc.ZPlane(z0=10.0, boundary_type='vacuum')
    bottom = openmc.ZPlane(z0=-10.0, boundary_type='vacuum')
    cell1 = openmc.Cell(fill=fuel, region=box1 & +bottom & -top)
    cell2 = openmc.Cell(fill=zr, region=~box1 & box2 & +bottom & -top)
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

    recti_mesh = openmc.RectilinearMesh()
    recti_mesh.x_grid = np.linspace(-7.5, 7.5, 18)
    recti_mesh.y_grid = np.linspace(-7.5, 7.5, 18)
    recti_mesh.z_grid = np.logspace(0, np.log10(7.5), 11)

    # Create filters
    reg_filters = [
        openmc.MeshFilter(mesh_1d),
        openmc.MeshFilter(mesh_2d),
        openmc.MeshFilter(mesh_3d),
        openmc.MeshFilter(recti_mesh)
    ]
    surf_filters = [
        openmc.MeshSurfaceFilter(mesh_1d),
        openmc.MeshSurfaceFilter(mesh_2d),
        openmc.MeshSurfaceFilter(mesh_3d),
        openmc.MeshSurfaceFilter(recti_mesh)
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
