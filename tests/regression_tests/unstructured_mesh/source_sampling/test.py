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

    dimen = 10
    size_hex = 20.0 / dimen

    cells = np.empty((dimen, dimen, dimen), dtype=object)
    surfaces = np.empty((dimen + 1, 3), dtype=object)

    geometry = openmc.Geometry()
    universe = openmc.Universe(universe_id=1, name="Contains all hexes")

    for i in range(0,dimen+1):
        coord = -dimen + i * size_hex
        surfaces[i][0] = openmc.XPlane(coord, name=f"X plane at {coord}")
        surfaces[i][1] = openmc.YPlane(coord, name=f"Y plane at {coord}")
        surfaces[i][2] = openmc.ZPlane(coord, name=f"Z plane at {coord}")

        surfaces[i][0].boundary_type = 'vacuum'
        surfaces[i][1].boundary_type = 'vacuum'
        surfaces[i][2].boundary_type = 'vacuum'

    for (k, j, i) in np.ndindex(cells.shape):
        cells[i][j][k] = openmc.Cell(name=("x = {}, y = {}, z = {}".format(i,j,k)))
        cells[i][j][k].region = +surfaces[i][0] & -surfaces[i+1][0] & \
                                +surfaces[j][1] & -surfaces[j+1][1] & \
                                +surfaces[k][2] & -surfaces[k+1][2]
        cells[i][j][k].fill = None

        universe.add_cell(cells[i][j][k])

    geometry = openmc.Geometry(universe)

    ### Settings ###
    settings = openmc.Settings()
    settings.run_mode = 'fixed source'
    settings.particles = 100
    settings.batches = 2

    mesh_filename = "test_mesh_tets.e"

    uscd_mesh = openmc.UnstructuredMesh(mesh_filename, 'libmesh')

    n_cells = len(geometry.get_all_cells())

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