import glob
from itertools import product
import os

import pytest
import numpy as np

import openmc
import openmc.lib

from tests import cdtemp
from tests.regression_tests import config
from subprocess import call

TETS_PER_VOXEL = 12


@pytest.fixture
def model():
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
    settings.particles = 5000
    settings.batches = 2

    return openmc.model.Model(geometry=geometry,
                              materials=materials,
                              settings=settings)


param_values = (['libmesh', 'moab'], # mesh libraries
                ['uniform', 'manual']) # Element weighting schemes

test_cases = []
for i, (lib, schemes) in enumerate(product(*param_values)):
    test_cases.append({'library' : lib,
                       'source_strengths' : schemes,
                       'inputs_true' : 'inputs_true{}.dat'.format(i)})

def ids(params):
    return f"{params['library']}-{params['source_strengths']}"

@pytest.mark.parametrize("test_cases", test_cases, ids=ids)
def test_unstructured_mesh_sampling(model, test_cases):
    # skip the test if the library is not enabled
    if test_cases['library'] == 'moab' and not openmc.lib._dagmc_enabled():
        pytest.skip("DAGMC (and MOAB) mesh not enabled in this build.")

    if test_cases['library'] == 'libmesh' and not openmc.lib._libmesh_enabled():
        pytest.skip("LibMesh is not enabled in this build.")

    # setup mesh source ###
    mesh_filename = "test_mesh_tets.e"

    uscd_mesh = openmc.UnstructuredMesh(mesh_filename, test_cases['library'])

    n_cells = len(model.geometry.get_all_cells())

    # set source weights according to test case
    if test_cases['source_strengths'] == 'uniform':
        vol_norm = True
        strengths = None
    elif test_cases['source_strengths'] == 'manual':
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
    model.settings.source = source

    with cdtemp(['test_mesh_tets.e', 'test_mesh_tets.exo']):
        model.export_to_xml()

        n_cells = len(model.geometry.get_all_cells())

        n_samples = 100000

        cell_counts = np.zeros(n_cells)

        # This model contains 1000 geometry cells. Each cell is a hex
        # corresponding to 12 of the tets. This test runs 10000 particles. This
        #  results in the following average for each cell
        average_in_hex = n_samples / n_cells

        openmc.lib.init([])

        sites = openmc.lib.sample_external_source(n_samples)
        cells = [openmc.lib.find_cell(s.r) for s in sites]

        openmc.lib.finalize()

        for c in cells:
            cell_counts[c[0]._index] += 1

        source_strengths = model.settings.source[0].space.strengths

        if source_strengths is not None:
            assert(cell_counts[0] > 0 and cell_counts[1] > 0)
            assert(cell_counts[0] > cell_counts[1])

            # counts for all other cells should be zero
            for i in range(2, len(cell_counts)):
                assert(cell_counts[i] == 0)

        else:
            # check that the average number of source sites in each cell
            # is within the expected deviation
            diff = np.abs(cell_counts - average_in_hex)

            assert(np.average(cell_counts) == average_in_hex) # this probably shouldn't be exact???
            assert((diff < 2*cell_counts.std()).sum() / diff.size >= 0.75)
            assert((diff < 6*cell_counts.std()).sum() / diff.size >= 0.97)
