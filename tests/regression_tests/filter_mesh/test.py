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

    # Initialize Meshes
    mesh_1d = openmc.RegularMesh()
    mesh_1d.dimension = [17]
    mesh_1d.lower_left = [-10.0]
    mesh_1d.upper_right = [10.0]

    mesh_2d = openmc.RegularMesh()
    mesh_2d.dimension = [17, 17]
    mesh_2d.lower_left = [-10.0, -10.0]
    mesh_2d.upper_right = [10.0, 10.0]

    mesh_3d = openmc.RegularMesh()
    mesh_3d.dimension = [17, 17, 17]
    mesh_3d.lower_left = [-10.0, -10.0, -183.00]
    mesh_3d.upper_right = [10.0, 10.0, 183.00]

    recti_mesh = openmc.RectilinearMesh()
    recti_mesh.x_grid = np.linspace(-10.0, 10.0, 18)
    recti_mesh.y_grid = np.linspace(-10.0, 10.0, 18)
    recti_mesh.z_grid = np.logspace(0, np.log10(183), 11)

    # Initialize the filters
    mesh_1d_filter = openmc.MeshFilter(mesh_1d)
    mesh_2d_filter = openmc.MeshFilter(mesh_2d)
    mesh_3d_filter = openmc.MeshFilter(mesh_3d)
    recti_mesh_filter = openmc.MeshFilter(recti_mesh)
    meshsurf_1d_filter = openmc.MeshSurfaceFilter(mesh_1d)
    meshsurf_2d_filter = openmc.MeshSurfaceFilter(mesh_2d)
    meshsurf_3d_filter = openmc.MeshSurfaceFilter(mesh_3d)
    recti_meshsurf_filter = openmc.MeshSurfaceFilter(recti_mesh)

    # Initialized the tallies
    tally = openmc.Tally(name='tally 1')
    tally.filters = [mesh_1d_filter]
    tally.scores = ['total']
    model.tallies.append(tally)

    tally = openmc.Tally(name='tally 2')
    tally.filters = [meshsurf_1d_filter]
    tally.scores = ['current']
    model.tallies.append(tally)

    tally = openmc.Tally(name='tally 3')
    tally.filters = [mesh_2d_filter]
    tally.scores = ['total']
    model.tallies.append(tally)

    tally = openmc.Tally(name='tally 4')
    tally.filters = [meshsurf_2d_filter]
    tally.scores = ['current']
    model.tallies.append(tally)

    tally = openmc.Tally(name='tally 5')
    tally.filters = [mesh_3d_filter]
    tally.scores = ['total']
    model.tallies.append(tally)

    tally = openmc.Tally(name='tally 6')
    tally.filters = [meshsurf_3d_filter]
    tally.scores = ['current']
    model.tallies.append(tally)

    tally = openmc.Tally(name='tally 7')
    tally.filters = [recti_mesh_filter]
    tally.scores = ['total']
    model.tallies.append(tally)

    tally = openmc.Tally(name='tally 8')
    tally.filters = [recti_meshsurf_filter]
    tally.scores = ['current']
    model.tallies.append(tally)

    return model


def test_filter_mesh(model):
    harness = HashedPyAPITestHarness('statepoint.5.h5', model)
    harness.main()
