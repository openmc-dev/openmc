from itertools import product

import pytest
import numpy as np

import openmc
import openmc.lib

from tests import cdtemp
from tests.testing_harness import PyAPITestHarness
from tests.regression_tests import config
from subprocess import call

TETS_PER_VOXEL = 12

# This test uses a geometry file that resembles a regular mesh.


@pytest.fixture
def model():
    openmc.reset_auto_ids()

    ### Materials ###
    materials = openmc.Materials()

    water_mat = openmc.Material(material_id=3, name="water")
    water_mat.add_nuclide("H1", 2.0)
    water_mat.add_nuclide("O16", 1.0)
    water_mat.set_density("atom/b-cm", 0.07416)
    materials.append(water_mat)

    ### Geometry ###
    # This test uses a geometry file that resembles a regular mesh.
    # 12 tets are used to match each voxel in the geometry.

    # create a regular mesh that matches the superimposed mesh
    regular_mesh = openmc.RegularMesh(mesh_id=10)
    regular_mesh.lower_left = (-10, -10, -10)
    regular_mesh.dimension = (10, 10, 10)
    regular_mesh.width = (2, 2, 2)

    root_cell, cells = regular_mesh.build_cells()

    geometry = openmc.Geometry(root=[root_cell])

    # set boundary conditions of the root cell
    for surface in root_cell.region.get_surfaces().values():
        surface.boundary_type = 'vacuum'

    ### Settings ###
    settings = openmc.Settings()
    settings.run_mode = 'fixed source'
    settings.particles = 100
    settings.batches = 2

    mesh_filename = "test_mesh_tets.e"
    uscd_mesh = openmc.UnstructuredMesh(mesh_filename, 'libmesh')

    # subtract one to account for root cell
    n_cells = len(geometry.get_all_cells()) - 1

    # set source weights according to test case
    vol_norm = False
    strengths = np.zeros(n_cells*TETS_PER_VOXEL)
    # set non-zero strengths only for the tets corresponding to the
    # first two geometric hex cells
    strengths[0:TETS_PER_VOXEL] = 10
    strengths[TETS_PER_VOXEL:2*TETS_PER_VOXEL] = 2

    # create the spatial distribution based on the mesh
    space = openmc.stats.MeshSpatial(uscd_mesh, strengths, vol_norm)

    energy = openmc.stats.Discrete(x=[15.e+06], p=[1.0])
    source = openmc.Source(space=space, energy=energy)
    settings.source = source
    return openmc.model.Model(geometry=geometry,
                              materials=materials,
                              settings=settings)

test_cases = []
for i, lib, in enumerate(['libmesh', 'moab']):
    test_cases.append({'library' : lib,
                       'inputs_true' : 'inputs_true{}.dat'.format(i)})

def ids(params):
    return params['library']

@pytest.mark.parametrize("test_cases", test_cases, ids=ids)
def test_unstructured_mesh_sampling(model, test_cases):

    model.settings.source[0].space.mesh.libaray = test_cases['library']

    harness = PyAPITestHarness('statepoint.2.h5', model, test_cases['inputs_true'])
    harness.main()