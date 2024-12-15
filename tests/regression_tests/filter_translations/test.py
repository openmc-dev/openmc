import numpy as np

import openmc
import pytest

from tests.testing_harness import PyAPITestHarness


@pytest.fixture
def model():

    model = openmc.model.Model()

    fuel = openmc.Material()
    fuel.set_density("g/cm3", 10.0)
    fuel.add_nuclide("U235", 1.0)
    zr = openmc.Material()
    zr.set_density("g/cm3", 1.0)
    zr.add_nuclide("Zr90", 1.0)
    model.materials.extend([fuel, zr])

    box1 = openmc.model.RectangularPrism(10.0, 10.0)
    box2 = openmc.model.RectangularPrism(20.0, 20.0, boundary_type="reflective")
    top = openmc.ZPlane(z0=10.0, boundary_type="vacuum")
    bottom = openmc.ZPlane(z0=-10.0, boundary_type="vacuum")
    cell1 = openmc.Cell(fill=fuel, region=-box1 & +bottom & -top)
    cell2 = openmc.Cell(fill=zr, region=+box1 & -box2 & +bottom & -top)
    model.geometry = openmc.Geometry([cell1, cell2])

    model.settings.batches = 5
    model.settings.inactive = 0
    model.settings.particles = 1000

    translation = np.array((10, -5, 0))

    llc = np.array([-9, -9, -9])
    urc = np.array([9, 9, 9])

    mesh_dims = (3, 4, 5)

    filters = []

    # un-translated meshes
    reg_mesh = openmc.RegularMesh()
    reg_mesh.dimension = mesh_dims
    reg_mesh.lower_left = llc
    reg_mesh.upper_right = urc

    filters.append(openmc.MeshFilter(reg_mesh))

    recti_mesh = openmc.RectilinearMesh()
    recti_mesh.x_grid = np.linspace(llc[0], urc[0], mesh_dims[0])
    recti_mesh.y_grid = np.linspace(llc[1], urc[1], mesh_dims[1])
    recti_mesh.z_grid = np.linspace(llc[2], urc[2], mesh_dims[2])

    filters.append(openmc.MeshFilter(recti_mesh))

    llc = np.array(llc - translation)
    urc = np.array(urc - translation)

    # translated meshes
    translated_reg_mesh = openmc.RegularMesh()
    translated_reg_mesh.dimension = mesh_dims
    translated_reg_mesh.lower_left = llc
    translated_reg_mesh.upper_right = urc

    filters.append(openmc.MeshFilter(translated_reg_mesh))
    filters[-1].translation = translation

    translated_recti_mesh = openmc.RectilinearMesh()
    translated_recti_mesh.x_grid = np.linspace(llc[0], urc[0], mesh_dims[0])
    translated_recti_mesh.y_grid = np.linspace(llc[1], urc[1], mesh_dims[1])
    translated_recti_mesh.z_grid = np.linspace(llc[2], urc[2], mesh_dims[2])

    filters.append(openmc.MeshFilter(translated_recti_mesh))
    filters[-1].translation = translation

    # Create tallies
    for f in filters:
        tally = openmc.Tally()
        tally.filters = [f]
        tally.scores = ["total"]
        model.tallies.append(tally)

    return model


def test_filter_mesh_translations(model):
    harness = PyAPITestHarness("statepoint.5.h5", model)
    harness.main()
