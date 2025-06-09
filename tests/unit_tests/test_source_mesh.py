from itertools import product
from pathlib import Path
from math import sqrt
import random

import pytest
import numpy as np
import openmc
import openmc.lib

from tests import cdtemp


###################
# MeshSpatial Tests
###################
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

    return openmc.Model(geometry=geometry,
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


###################
# MeshSource tests
###################
@pytest.fixture
def void_model():
    """
    A void model containing a single box
    """
    model = openmc.Model()

    box = openmc.model.RectangularParallelepiped(*[-10, 10]*3, boundary_type='vacuum')
    model.geometry = openmc.Geometry([openmc.Cell(region=-box)])

    model.settings.particles = 100
    model.settings.batches = 10
    model.settings.run_mode = 'fixed source'

    return model


@pytest.mark.parametrize('mesh_type', ('rectangular', 'cylindrical'))
def test_mesh_source_independent(run_in_tmpdir, void_model, mesh_type):
    """
    A void model containing a single box
    """
    model = void_model

    # define a 2 x 2 x 2 mesh
    if mesh_type == 'rectangular':
        mesh = openmc.RegularMesh.from_domain(model.geometry, (2, 2, 2))
    elif mesh_type == 'cylindrical':
        mesh = openmc.CylindricalMesh.from_domain(model.geometry, (1, 4, 2))

    energy = openmc.stats.Discrete([1.e6], [1.0])

    # create sources with only one non-zero strength for the source in the mesh
    # voxel occupying the lowest octant. Direct source particles straight out of
    # the problem from there. This demonstrates that
    # 1) particles are only being sourced within the intented mesh voxel based
    #    on source strength
    # 2) particles are respecting the angle distributions assigned to each voxel
    sources = np.empty(mesh.dimension, dtype=openmc.SourceBase)
    centroids = mesh.centroids
    x, y, z = np.swapaxes(mesh.centroids, -1, 0)
    for i, j, k in mesh.indices:
        # mesh.indices is currently one-indexed, adjust for Python arrays
        ijk = (i-1, j-1, k-1)

        # get the centroid of the ijk mesh element and use it to set the
        # direction of the source directly out of the problem
        centroid = centroids[ijk]
        vec = np.sign(centroid, dtype=float)
        vec /= np.linalg.norm(vec)
        angle = openmc.stats.Monodirectional(vec)
        sources[ijk] = openmc.IndependentSource(energy=energy, angle=angle, strength=0.0)

    # create and apply the mesh source
    mesh_source = openmc.MeshSource(mesh, sources)
    model.settings.source = mesh_source

    # tally the flux on the mesh
    mesh_filter = openmc.MeshFilter(mesh)
    tally = openmc.Tally()
    tally.filters = [mesh_filter]
    tally.scores = ['flux']

    model.tallies = openmc.Tallies([tally])

    # for each element, set a single-non zero source with particles
    # traveling out of the mesh (and geometry) w/o crossing any other
    # mesh elements
    for flat_index, (i, j, k) in enumerate(mesh.indices):
        ijk = (i-1, j-1, k-1)
        # zero-out all source strengths and set the strength
        # on the element of interest
        mesh_source.strength = 0.0
        mesh_source.sources[flat_index].strength = 1.0

        sp_file = model.run()

        with openmc.StatePoint(sp_file) as sp:
            tally_out = sp.get_tally(id=tally.id)
            mean = tally_out.get_reshaped_data(expand_dims=True)

        # remove nuclides and scores axes
        mean = mean[..., 0, 0]
        # the mesh elment with a non-zero source strength should have a value
        assert mean[ijk] != 0
        # all other values should be zero
        mean[ijk] = 0
        assert np.all(mean == 0), f'Failed on index {ijk} with centroid {mesh.centroids[ijk]}'

        # test roundtrip
        xml_model = openmc.Model.from_model_xml()
        xml_source = xml_model.settings.source[0]
        assert isinstance(xml_source, openmc.MeshSource)
        assert xml_source.strength == 1.0
        assert isinstance(xml_source.mesh, type(mesh_source.mesh))
        assert xml_source.mesh.dimension == mesh_source.mesh.dimension
        assert xml_source.mesh.id == mesh_source.mesh.id
        assert len(xml_source.sources) == len(mesh_source.sources)

    # check strength adjustment methods
    assert mesh_source.strength == 1.0
    mesh_source.strength = 100.0
    assert mesh_source.strength == 100.0

    mesh_source.normalize_source_strengths()
    assert mesh_source.strength == 1.0


@pytest.mark.parametrize("library", ('moab', 'libmesh'))
def test_umesh_source_independent(run_in_tmpdir, request, void_model, library):
    import openmc.lib
    # skip the test if the library is not enabled
    if library == 'moab' and not openmc.lib._dagmc_enabled():
        pytest.skip("DAGMC (and MOAB) mesh not enabled in this build.")

    if library == 'libmesh' and not openmc.lib._libmesh_enabled():
        pytest.skip("LibMesh is not enabled in this build.")

    model = void_model

    mesh_filename = Path(request.fspath).parent / "test_mesh_tets.e"
    uscd_mesh = openmc.UnstructuredMesh(mesh_filename, library)
    ind_source = openmc.IndependentSource()
    n_elements = 12_000
    model.settings.source = openmc.MeshSource(uscd_mesh, n_elements*[ind_source])
    model.export_to_model_xml()
    with openmc.lib.run_in_memory():
        openmc.lib.simulation_init()
        sites = openmc.lib.sample_external_source(10)
        openmc.lib.statepoint_write('statepoint.h5')

    with openmc.StatePoint('statepoint.h5') as sp:
        uscd_mesh = sp.meshes[uscd_mesh.id]

    # ensure at least that all sites are inside the mesh
    bounding_box = uscd_mesh.bounding_box
    for site in sites:
        assert site.r in bounding_box


def test_mesh_source_constraints(run_in_tmpdir):
    """Test application of constraints to underlying mesh element sources"""

    # Create simple model with two cells
    m1 = openmc.Material()
    m1.add_nuclide('H1', 1.0)
    m2 = m1.clone()
    sph = openmc.Sphere(r=100, boundary_type='vacuum')
    box1 = openmc.model.RectangularParallelepiped(-1, 0, -1, 1, -1, 1)
    box2 = openmc.model.RectangularParallelepiped(0, 2, -1, 1, -1, 1)
    cell1 = openmc.Cell(fill=m1, region=-box1)
    cell2 = openmc.Cell(fill=m2, region=-box2)
    outer = openmc.Cell(region=-sph & (+box1 | +box2))
    model = openmc.Model()
    model.geometry = openmc.Geometry([cell1, cell2, outer])

    # Define a mesh covering the two cells: the first mesh element contains
    # cell1 (-1 < x < 0) and the second element contains cells2 (0 < x < 2)
    mesh = openmc.RegularMesh()
    mesh.lower_left = (-3., -1., -1.)
    mesh.upper_right = (3., 1., 1.)
    mesh.dimension = (2, 1, 1)

    # Define a mesh source with a randomly chosen probability
    p = random.random()
    src1 = openmc.IndependentSource(strength=p, constraints={'domains': [cell1]})
    src2 = openmc.IndependentSource(strength=1 - p, constraints={'domains': [cell2]})
    model.settings.source = openmc.MeshSource(mesh, [src1, src2])

    # Finish settings and export
    model.settings.particles = 100
    model.settings.batches = 1
    model.export_to_model_xml()

    with openmc.lib.run_in_memory():
        # Sample sites from the source
        sites = openmc.lib.sample_external_source(N := 1000)

        # Check that all sites are either in cell1 or cell2
        xs = np.array([s.r[0] for s in sites])
        assert (xs >= -1.0).all()
        assert (xs <= 2.0).all()

        # Check that the correct percentage of the sites are in cell1
        sigma = sqrt(p*(1- p)/N)
        frac = xs[(-1.0 <= xs) & (xs <= 0.0)].size / N
        assert frac == pytest.approx(p, abs=5*sigma)
