from itertools import product
from pathlib import Path
from subprocess import call

import pytest
import numpy as np
import openmc
import openmc.lib

from tests import cdtemp
from tests.regression_tests import config


TETS_PER_VOXEL = 12

# This test uses a geometry file with cells that match a regular mesh. Each cell
# in the geometry corresponds to 12 tetrahedra in the unstructured mesh file.
@pytest.fixture
def model():
    openmc.reset_auto_ids()

    ### Materials ###
    materials = openmc.Materials()

    water_mat = openmc.Material(name="water")
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

    root_cell, _ = regular_mesh.build_cells(bc=['vacuum']*6)

    geometry = openmc.Geometry(root=[root_cell])

    ### Settings ###
    settings = openmc.Settings()
    settings.run_mode = 'fixed source'
    settings.particles = 100
    settings.batches = 2

    return openmc.model.Model(geometry=geometry,
                              materials=materials,
                              settings=settings)

### Setup test cases ###
param_values = (['libmesh', 'moab'], # mesh libraries
                ['uniform', 'manual']) # Element weighting schemes

test_cases = []
for i, (lib, schemes) in enumerate(product(*param_values)):
    test_cases.append({'library' : lib,
                       'source_strengths' : schemes})

def ids(params):
    """Test naming function for clarity"""
    return f"{params['library']}-{params['source_strengths']}"

@pytest.mark.parametrize("test_cases", test_cases, ids=ids)
def test_unstructured_mesh_sampling(model, request, test_cases):
    # skip the test if the library is not enabled
    if test_cases['library'] == 'moab' and not openmc.lib._dagmc_enabled():
        pytest.skip("DAGMC (and MOAB) mesh not enabled in this build.")

    if test_cases['library'] == 'libmesh' and not openmc.lib._libmesh_enabled():
        pytest.skip("LibMesh is not enabled in this build.")

    # setup mesh source ###
    mesh_filename = Path(request.fspath).parent / "test_mesh_tets.e"
    uscd_mesh = openmc.UnstructuredMesh(mesh_filename, test_cases['library'])

    # subtract one to account for root cell produced by RegularMesh.build_cells
    n_cells = len(model.geometry.get_all_cells()) - 1

    # set source weights according to test case
    if test_cases['source_strengths'] == 'uniform':
        vol_norm = True
        strengths = None
    elif test_cases['source_strengths'] == 'manual':
        vol_norm = False
        # assign random weights
        strengths = np.random.rand(n_cells*TETS_PER_VOXEL)

    # create the spatial distribution based on the mesh
    space = openmc.stats.MeshSpatial(uscd_mesh, strengths, vol_norm)

    energy = openmc.stats.Discrete(x=[15.e+06], p=[1.0])
    source = openmc.IndependentSource(space=space, energy=energy)
    model.settings.source = source

    with cdtemp([mesh_filename]):
        model.export_to_xml()

        n_measurements = 100
        n_samples = 1000

        cell_counts = np.zeros((n_cells, n_measurements))

        # This model contains 1000 geometry cells. Each cell is a hex
        # corresponding to 12 of the tets. This test runs 1000 samples. This
        #  results in the following average for each cell
        openmc.lib.init([])

        # perform many sets of samples and track counts for each cell
        for m in range(n_measurements):
            sites = openmc.lib.sample_external_source(n_samples)
            cells = [openmc.lib.find_cell(s.r) for s in sites]

            for c in cells:
                # subtract one from index to account for root cell
                cell_counts[c[0]._index - 1, m] += 1

        # make sure particle transport is successful
        openmc.lib.run()
        openmc.lib.finalize()

        # normalize cell counts to get sampling frequency per particle
        cell_counts /= n_samples

        # get the mean and std. dev. of the cell counts
        mean = cell_counts.mean(axis=1)
        std_dev = cell_counts.std(axis=1)

        if test_cases['source_strengths'] == 'uniform':
            exp_vals = np.ones(n_cells) / n_cells
        else:
            # sum up the source strengths for each tet, these are the expected true mean
            # of the sampling frequency for that cell
            exp_vals = strengths.reshape(-1, 12).sum(axis=1) / sum(strengths)

        diff = np.abs(mean - exp_vals)
        assert((diff < 2*std_dev).sum() / diff.size >= 0.95)
        assert((diff < 6*std_dev).sum() / diff.size >= 0.997)


def test_strengths_size_failure(request, model):
    # setup mesh source ###
    mesh_filename = Path(request.fspath).parent / "test_mesh_tets.e"
    uscd_mesh = openmc.UnstructuredMesh(mesh_filename, 'libmesh')

    # intentionally incorrectly sized to trigger an error
    n_cells = len(model.geometry.get_all_cells())
    strengths = np.random.rand(n_cells*TETS_PER_VOXEL)

    # create the spatial distribution based on the mesh
    space = openmc.stats.MeshSpatial(uscd_mesh, strengths)

    energy = openmc.stats.Discrete(x=[15.e+06], p=[1.0])
    source = openmc.IndependentSource(space=space, energy=energy)
    model.settings.source = source

    # skip the test if unstructured mesh is not available
    if not openmc.lib._libmesh_enabled():
        if openmc.lib._dagmc_enabled():
            source.space.mesh.library = 'moab'
        else:
            pytest.skip("Unstructured mesh support unavailable.")

    # make sure that an incorrrectly sized strengths array causes a failure
    source.space.strengths = source.space.strengths[:-1]

    mesh_filename = Path(request.fspath).parent / source.space.mesh.filename

    with pytest.raises(RuntimeError, match=r'strengths array'), cdtemp([mesh_filename]):
        model.export_to_xml()
        openmc.run()

def test_roundtrip(run_in_tmpdir, model, request):
    if not openmc.lib._libmesh_enabled() and not openmc.lib._dagmc_enabled():
        pytest.skip("Unstructured mesh is not enabled in this build.")

    mesh_filename = Path(request.fspath).parent / 'test_mesh_tets.e'
    ucd_mesh = openmc.UnstructuredMesh(mesh_filename, library='libmesh')

    if not openmc.lib._libmesh_enabled():
        ucd_mesh.library = 'moab'

    n_cells = len(model.geometry.get_all_cells())

    space_out = openmc.MeshSpatial(ucd_mesh)
    space_out.strengths = np.random.rand(n_cells*TETS_PER_VOXEL)
    model.settings.source = openmc.IndependentSource(space=space_out)

    # write out the model
    model.export_to_xml()

    model_in = openmc.Model.from_xml()

    space_in = model_in.settings.source[0].space

    np.testing.assert_equal(space_out.strengths, space_in.strengths)

    assert space_in.mesh.id == space_out.mesh.id
    assert space_in.volume_normalized == space_out.volume_normalized
