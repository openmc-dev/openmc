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

class UnstructuredMeshSourceTest(PyAPITestHarness):
    def __init__(self, statepoint_name, model, inputs_true, schemes):
        super().__init__(statepoint_name, model, inputs_true)
        self.schemes = schemes

    def _run_openmc(self):
        kwargs = {'openmc_exec' : config['exe'],
                  'event_based' : config['event'],
                  'tracks' : "True"}

        if config['mpi']:
            kwargs['mpi_args'] = [config['mpiexec'], '-n', config['mpi_np']]

        openmc.run(**kwargs)

    def _compare_results(self):
        # There are 10000 particles and 1000 hexes, this leads to the average
        # shown below
        average_in_hex = 10.0

        # Load in tracks
        if config['mpi']:
            openmc.Tracks.combine(glob.glob('tracks_p*.h5'))
        
        tracks = openmc.Tracks(filepath='tracks.h5')
        tracks_born = np.empty((len(tracks), 1))

        instances = np.zeros(1000)

        # loop over the tracks and get data
        for i in range(0, len(tracks)):
            tracks_born[i] = tracks[i].particle_tracks[0].states['cell_id'][0]
            instances[int(tracks_born[i])-1] += 1

        if self.schemes == "file":
            assert(instances[0] > 0 and instances[1] > 0)
            assert(instances[0] > instances[1])

            for i in range(2, len(instances)):
                assert(instances[i] == 0)

        else:
            assert(np.average(instances) == average_in_hex)
            assert(np.std(instances) < np.average(instances))
            assert(np.amax(instances) < np.average(instances)+6*np.std(instances))


    def _cleanup(self):
        super()._cleanup()
        output = glob.glob('track*.h5')
        output += glob.glob('tally*.e')
        for f in output:
            if os.path.exists(f):
                os.remove(f)

param_values = (['libmesh', 'moab'], # mesh libraries
                ['volume', 'file']) # Element weighting schemes

test_cases = []
for i, (lib, schemes) in enumerate(product(*param_values)):
    test_cases.append({'library' : lib,
                       'schemes' : schemes,
                       'inputs_true' : 'inputs_true{}.dat'.format(i)})

@pytest.mark.parametrize("test_opts", test_cases)
def test_unstructured_mesh(test_opts):

    openmc.reset_auto_ids()

    # skip the test if the library is not enabled
    if test_opts['library'] == 'moab' and not openmc.lib._dagmc_enabled():
        pytest.skip("DAGMC (and MOAB) mesh not enabled in this build.")

    if test_opts['library'] == 'libmesh' and not openmc.lib._libmesh_enabled():
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

    mesh_filename = "test_mesh_tets.e"

    uscd_mesh = openmc.UnstructuredMesh(mesh_filename, test_opts['library'])

    ### Tallies ###

    # create tallies
    tallies = openmc.Tallies()

    tally1 = openmc.Tally(1)
    tally1.scores = ['scatter', 'total', 'absorption']
    # Export tallies
    tallies = openmc.Tallies([tally1])
    tallies.export_to_xml()

    ### Settings ###
    settings = openmc.Settings()
    settings.run_mode = 'fixed source'
    settings.particles = 5000
    settings.batches = 2

    settings.max_tracks = 10000

    # source setup
    if test_opts['schemes'] == 'volume':
        space = openmc.stats.MeshSpatial(volume_normalized=True, mesh=uscd_mesh)
    elif test_opts['schemes'] == 'file':
        array = np.zeros(12000)
        for i in range(0, 12):
            array[i] = 10
            array[i+12] = 2
        space = openmc.stats.MeshSpatial(volume_normalized=False, strengths=array, mesh=uscd_mesh)

    energy = openmc.stats.Discrete(x=[15.e+06], p=[1.0])
    source = openmc.Source(space=space, energy=energy)
    settings.source = source

    model = openmc.model.Model(geometry=geometry,
                               materials=materials,
                               tallies=tallies,
                               settings=settings)

    harness = UnstructuredMeshSourceTest('statepoint.2.h5',
                                         model,
                                         test_opts['inputs_true'],
                                         test_opts['schemes'])
    harness.main()
