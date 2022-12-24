import glob
from itertools import product
import os

import pytest
import numpy as np

import openmc
import openmc.lib

from tests.testing_harness import PyAPITestHarness
from tests.regression_tests import config
from subprocess import call

TETS_PER_VOXEL = 12

# This test uses a geometry file that resembles a regular mesh.
# 12 tets are used to match each voxel in the geometry.

class UnstructuredMeshSourceTest(PyAPITestHarness):
    def __init__(self, statepoint_name, model, inputs_true):
        super().__init__(statepoint_name, model, inputs_true)

    def _run_openmc(self):
        kwargs = {'openmc_exec' : config['exe'],
                  'event_based' : config['event'],
                  'tracks' : "True"}

        if config['mpi']:
            kwargs['mpi_args'] = [config['mpiexec'], '-n', config['mpi_np']]

        openmc.run(**kwargs)

    def _compare_results(self):
        # This model contains 1000 geometry cells. Each cell is a hex
        # corresponding to 12 of the tets. This test runs 10000 particles. This
        #  results in the following average for each cell

        # we can compute this based on the number of particles run in the simulation
        average_in_hex = 10.0

        # Load in tracks
        if config['mpi']:
            openmc.Tracks.combine(glob.glob('tracks_p*.h5'))

        tracks = openmc.Tracks(filepath='tracks.h5')
        tracks_born = np.empty((len(tracks), 1))

        # create an array with an entry for each geometric cell
        cell_counts = np.zeros(1000)

        # loop over the tracks and get data
        for i in range(0, len(tracks)):
            # get the initial cell ID of the track, and assign it for the tracks_born array
            tracks_born[i] = tracks[i].particle_tracks[0].states['cell_id'][0]
            # increment the cell_counts entry for this cell_id
            cell_counts[int(tracks_born[i])-1] += 1

        source_strengths = self._model.settings.source[0].space.strengths

        if source_strengths is not None:
            assert(cell_counts[0] > 0 and cell_counts[1] > 0)
            assert(cell_counts[0] > cell_counts[1])

            # counts for all other cells should be zero
            for i in range(2, len(cell_counts)):
                assert(cell_counts[i] == 0)

        else:
            # check that the average number of source sites in each cell
            diff = np.abs(cell_counts - average_in_hex)

            assert(np.average(cell_counts) == average_in_hex) # this probably shouldn't be exact???
            assert((diff < 2*cell_counts.std()).sum() / diff.size >= 0.75)
            assert((diff < 6*cell_counts.std()).sum() / diff.size >= 0.97)

    def _cleanup(self):
        super()._cleanup()
        output = glob.glob('track*.h5')
        output += glob.glob('tally*.e')
        for f in output:
            if os.path.exists(f):
                os.remove(f)

param_values = (['libmesh', 'moab'], # mesh libraries
                ['uniform', 'manual']) # Element weighting schemes

test_cases = []
for i, (lib, schemes) in enumerate(product(*param_values)):
    test_cases.append({'library' : lib,
                       'source_strengths' : schemes,
                       'inputs_true' : 'inputs_true{}.dat'.format(i)})

@pytest.mark.parametrize("test_cases", test_cases)
def test_unstructured_mesh_sampling(test_cases):
    openmc.reset_auto_ids()

    # skip the test if the library is not enabled
    if test_cases['library'] == 'moab' and not openmc.lib._dagmc_enabled():
        pytest.skip("DAGMC (and MOAB) mesh not enabled in this build.")

    if test_cases['library'] == 'libmesh' and not openmc.lib._libmesh_enabled():
        pytest.skip("LibMesh is not enabled in this build.")

    ### Materials ###
    materials = openmc.Materials()

    water_mat = openmc.Material(material_id=3, name="water")
    water_mat.add_nuclide("H1", 2.0)
    water_mat.add_nuclide("O16", 1.0)
    water_mat.set_density("atom/b-cm", 0.07416)
    materials.append(water_mat)

    materials.export_to_xml()

    ### Geometry ###
    dimen = 10
    size_hex = 20.0 / dimen

    ### Geometry ###
    cell = np.empty((dimen, dimen, dimen), dtype=object)
    surfaces = np.empty((dimen + 1, 3), dtype=object)

    geometry = openmc.Geometry()
    universe = openmc.Universe(universe_id=1, name="Contains all hexes")

    for i in range(0,dimen+1):
        coord = -10.0 + i * size_hex
        surfaces[i][0] = openmc.XPlane(coord, name=f"X plane at {coord}")
        surfaces[i][1] = openmc.YPlane(coord, name=f"Y plane at {coord}")
        surfaces[i][2] = openmc.ZPlane(coord, name=f"Z plane at {coord}")

        surfaces[i][0].boundary_type = 'vacuum'
        surfaces[i][1].boundary_type = 'vacuum'
        surfaces[i][2].boundary_type = 'vacuum'

    for k in range(0,dimen):
        for j in range(0,dimen):
            for i in range(0,dimen):
                cell[i][j][k] = openmc.Cell(name=("x = {}, y = {}, z = {}".format(i,j,k)))
                cell[i][j][k].region = +surfaces[i][0] & -surfaces[i+1][0] & \
                                       +surfaces[j][1] & -surfaces[j+1][1] & \
                                       +surfaces[k][2] & -surfaces[k+1][2]
                cell[i][j][k].fill = None
                universe.add_cell(cell[i][j][k])

    geometry = openmc.Geometry(universe)

    ### Settings ###
    settings = openmc.Settings()
    settings.run_mode = 'fixed source'
    settings.particles = 5000
    settings.batches = 2

    settings.max_tracks = settings.particles * settings.batches

    ### Source ###
    mesh_filename = "test_mesh_tets.e"

    uscd_mesh = openmc.UnstructuredMesh(mesh_filename, test_cases['library'])

    # set source weights according to test case
    if test_cases['source_strengths'] == 'uniform':
        vol_norm = True
        strengths = None

    elif test_cases['source_strengths'] == 'manual':
        vol_norm = False
        strengths = np.zeros(12000)
        # set non-zero strengths only for the tets corresponding to the
        # first two geometric hex cells
        strengths[0:12] = 10
        strengths[12:24] = 2

    # create the spatial distribution based on the mesh
    space = openmc.stats.MeshSpatial(uscd_mesh, strengths, vol_norm)

    energy = openmc.stats.Discrete(x=[15.e+06], p=[1.0])
    source = openmc.Source(space=space, energy=energy)
    settings.source = source

    model = openmc.model.Model(geometry=geometry,
                               materials=materials,
                               settings=settings)
    harness = UnstructuredMeshSourceTest('statepoint.2.h5',
                                         model,
                                         test_cases['inputs_true'])
    harness.main()