from pathlib import Path

import openmc
import pytest

from tests.testing_harness import TestHarness
from tests.regression_tests import config
from tests import cdtemp

class StatepointRestartTestHarness(TestHarness):
    def __init__(self, final_sp, restart_sp):
        super().__init__(final_sp)
        self._restart_sp = restart_sp

    def execute_test(self):
        """Run OpenMC with the appropriate arguments and check the outputs."""
        try:
            self._run_openmc()
            self._test_output_created()
            results = self._get_results()
            self._write_results(results)
            self._compare_results()

            self._run_openmc_restart()
            self._test_output_created()
            results = self._get_results()
            self._write_results(results)
            self._compare_results()
        finally:
            self._cleanup()

    def update_results(self):
        """Update the results_true using the current version of OpenMC."""
        try:
            self._run_openmc()
            self._test_output_created()
            results = self._get_results()
            self._write_results(results)
            self._overwrite_results()
        finally:
            self._cleanup()

    def _run_openmc_restart(self):
        # Get the name of the statepoint file.
        statepoint = list(Path.cwd().glob(self._restart_sp))
        assert len(statepoint) == 1
        statepoint = statepoint[0]

        # Run OpenMC
        if config['mpi']:
            mpi_args = [config['mpiexec'], '-n', config['mpi_np']]
            openmc.run(restart_file=statepoint, openmc_exec=config['exe'],
                       mpi_args=mpi_args)
        else:
            openmc.run(openmc_exec=config['exe'], restart_file=statepoint)


def test_statepoint_restart():
    harness = StatepointRestartTestHarness('statepoint.10.h5',
                                           'statepoint.07.h5')
    harness.main()


def test_batch_check(request):
    xmls = list(request.path.parent.glob('*.xml'))

    with cdtemp(xmls):
        model = openmc.Model.from_xml()
        model.settings.particles = 100
        # run the model
        sp_file = model.run()

        # run a restart with the resulting statepoint
        # and the settings unchanged
        with pytest.raises(RuntimeError, match='is smaller than the number of batches'):
            model.run(restart_file=sp_file)

        # update the number of batches and run again
        model.settings.batches = 15
        model.settings.statepoint = {}
        sp_file = model.run(restart_file=sp_file)

        sp = openmc.StatePoint(sp_file)
        assert sp.n_batches == 15
