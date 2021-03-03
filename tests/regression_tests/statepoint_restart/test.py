import glob
import os

import openmc

from tests.testing_harness import TestHarness
from tests.regression_tests import config


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
        statepoint = glob.glob(os.path.join(os.getcwd(), self._restart_sp))
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
