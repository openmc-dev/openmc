import glob
import os

import openmc
import pytest

from tests.testing_harness import PyAPITestHarness
from tests.regression_tests import config


@pytest.fixture
def model():

    # Materials
    mat = openmc.Material()
    mat.set_density('g/cm3', 4.5)
    mat.add_nuclide('U235', 1.0)
    materials = openmc.Materials([mat])

    # Geometry
    sph = openmc.Sphere(r=10.0, boundary_type='vacuum')
    cell = openmc.Cell(fill=mat, region=-sph)
    geometry = openmc.Geometry([cell])

    # Settings
    settings = openmc.Settings()
    settings.run_mode = 'eigenvalue'
    settings.batches = 10
    settings.inactive = 5
    settings.particles = 200
    # Choose a sufficiently low threshold to trigger after more than 10 batches.
    # 0.004 seems to take 13 batches.
    settings.keff_trigger = {'type': 'std_dev', 'threshold': 0.004}
    settings.trigger_max_batches = 1000
    settings.trigger_batch_interval = 1
    settings.trigger_active = True
    settings.verbosity = 1 # to test that this works even with no output
 
    # Tallies
    t = openmc.Tally()
    t.scores = ['flux']
    tallies = openmc.Tallies([t])
 
    # Put it all together
    model = openmc.model.Model(materials=materials,
                               geometry=geometry,
                               settings=settings,
                               tallies=tallies)
    return model


class TriggerStatepointRestartTestHarness(PyAPITestHarness):
    def __init__(self, statepoint, model=None):
        super().__init__(statepoint, model)
        self._restart_sp = None
        self._final_sp = None
        # store the statepoint filename pattern separately to sp_name so we can reuse it
        self._sp_pattern = self._sp_name

    def _test_output_created(self):
        """Make sure statepoint files have been created."""
        spfiles = sorted(glob.glob(self._sp_pattern))
        assert len(spfiles) == 2, \
            'Two statepoint files should have been created'
        if not self._final_sp:
            # First non-restart run
            self._restart_sp = spfiles[0]
            self._final_sp = spfiles[1]
        else:
            # Second restart run
            assert spfiles[1] == self._final_sp, \
                'Final statepoint names were different'
        # Use the final_sp as the sp_name for the 'standard' results tests
        self._sp_name = self._final_sp

    def execute_test(self):
        """
        Perform initial and restart runs using the model.run method,
        Check all inputs and outputs which should be the same as those
        generated using the normal PyAPITestHarness update methods.
        """
        try:
            args = {'openmc_exec': config['exe'], 'event_based': config['event']}
            if config['mpi']:
                args['mpi_args'] = [config['mpiexec'], '-n', config['mpi_np']]

            # First non-restart run
            spfile = self._model.run(**args)
            sp_batchno_1 = 0
            print('Last sp file: %s' % spfile)
            assert spfile
            with openmc.StatePoint(spfile) as sp:
                 sp_batchno_1 = sp.current_batch
                 k_combined_1 = sp.k_combined
            assert sp_batchno_1 > 10
            print('Last batch no = %d' % sp_batchno_1)
            self._write_inputs(self._get_inputs())
            self._compare_inputs()
            self._test_output_created()
            self._write_results(self._get_results())
            self._compare_results()

            # Second restart run
            restart_spfile = glob.glob(os.path.join(os.getcwd(), self._restart_sp))
            assert len(restart_spfile) == 1
            args['restart_file'] = restart_spfile[0]
            spfile = self._model.run(**args)
            sp_batchno_2 = 0
            assert spfile
            with openmc.StatePoint(spfile) as sp:
                 sp_batchno_2 = sp.current_batch
                 k_combined_2 = sp.k_combined
            assert sp_batchno_2 > 10
            assert sp_batchno_1 == sp_batchno_2, \
                'Different final batch number after restart'
            # need str() here as uncertainties.ufloat instances are always different
            assert str(k_combined_1) == str(k_combined_2), \
                'Different final k_combined after restart'
            self._write_inputs(self._get_inputs())
            self._compare_inputs()
            self._test_output_created()
            self._write_results(self._get_results())
            self._compare_results()
        finally:
            self._cleanup()


def test_trigger_statepoint_restart(model):
    # Assuming we converge within 1000 batches, the statepoint filename
    # should include the batch number padded by at least one '0'.
    harness = TriggerStatepointRestartTestHarness('statepoint.0*.h5', model)
    harness.main()
