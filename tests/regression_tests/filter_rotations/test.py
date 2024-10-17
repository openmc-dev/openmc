import numpy as np

import openmc
import pytest

from tests.testing_harness import PyAPITestHarness


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

    rotation = np.array((180, 0, 0))

    llc = np.array([-9, -9, -9])
    urc = np.array([9, 9, 9])

    mesh_dims = (3, 4, 5)

    filters = []

    # un-rotated meshes
    reg_mesh = openmc.RegularMesh()
    reg_mesh.dimension = mesh_dims
    reg_mesh.lower_left = llc
    reg_mesh.upper_right = urc

    filters.append(openmc.MeshFilter(reg_mesh))


    # rotated meshes
    rotated_reg_mesh = openmc.RegularMesh()
    rotated_reg_mesh.dimension = mesh_dims
    rotated_reg_mesh.lower_left = llc
    rotated_reg_mesh.upper_right = urc

    filters.append(openmc.MeshFilter(rotated_reg_mesh))
    filters[-1].rotation = rotation

    # Create tallies
    for f in filters:
        tally = openmc.Tally()
        tally.filters = [f]
        tally.scores = ['total']
        model.tallies.append(tally)

    return model


def test_filter_mesh_rotations(model):
    harness = PyAPITestHarness('statepoint.5.h5', model)
    harness.main()
