from itertools import product

from pathlib import Path
import pytest
import numpy as np

import openmc
import openmc.lib

from tests import cdtemp
from tests.testing_harness import PyAPITestHarness

TETS_PER_VOXEL = 12

# This test uses a geometry file with cells that match a regular mesh. Each cell
# in the geometry corresponds to 12 tetrahedra in the unstructured mesh file.
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
    # create a regular mesh that matches the superimposed mesh
    regular_mesh = openmc.RegularMesh(mesh_id=10)
    regular_mesh.lower_left = (-10, -10, -10)
    regular_mesh.dimension = (10, 10, 10)
    regular_mesh.width = (2, 2, 2)

    root_cell, _ = regular_mesh.build_cells()

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

    # set source weights manually so the C++ that checks the
    # size of the array is executed
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

@pytest.mark.parametrize("test_cases", test_cases, ids=lambda p: p['library'])
def test_unstructured_mesh_sampling(model, test_cases):
    # skip the test if the library is not enabled
    if test_cases['library'] == 'moab' and not openmc.lib._dagmc_enabled():
        pytest.skip("DAGMC (and MOAB) mesh not enabled in this build.")

    if test_cases['library'] == 'libmesh' and not openmc.lib._libmesh_enabled():
        pytest.skip("LibMesh is not enabled in this build.")

    model.settings.source[0].space.mesh.libaray = test_cases['library']

    harness = PyAPITestHarness('statepoint.2.h5', model, test_cases['inputs_true'])
    harness.main()


def test_strengths_size_failure(request, model):
    mesh_source = model.settings.source[0]

    # skip the test if unstructured mesh is not available
    if not openmc.lib._libmesh_enabled():
        if openmc.lib._dagmc_enabled():
            mesh_source.space.mesh.library = 'moab'
        else:
            pytest.skip("Unstructured mesh support unavailable.")

    # make sure that an incorrrectly sized strengths array causes a failure
    mesh_source.space.strengths = mesh_source.space.strengths[:-1]

    mesh_filename = Path(request.fspath).parent / mesh_source.space.mesh.filename

    with pytest.raises(RuntimeError, match=r'strengths array'), cdtemp([mesh_filename]):
        model.export_to_xml()
        openmc.run()
