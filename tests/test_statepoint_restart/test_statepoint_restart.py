#!/usr/bin/env python

import glob
import os
import sys
sys.path.insert(0, os.pardir)
from testing_harness import TestHarness
from openmc.statepoint import StatePoint
from openmc.executor import Executor


class StatepointRestartTestHarness(TestHarness):
    def __init__(self, final_sp, restart_sp, tallies_present=False):
        super(StatepointRestartTestHarness, self).__init__(final_sp,
                                                           tallies_present)
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
        executor = Executor()

        if self._opts.mpi_exec is not None:
            returncode = executor.run_simulation(mpi_procs=self._opts.mpi_np,
                                                 restart_file=statepoint,
                                                 openmc_exec=self._opts.exe,
                                                 mpi_exec=self._opts.mpi_exec)

        else:
            returncode = executor.run_simulation(openmc_exec=self._opts.exe,
                                                 restart_file=statepoint)

        assert returncode == 0, 'OpenMC did not exit successfully.'


if __name__ == '__main__':
    harness = StatepointRestartTestHarness('statepoint.10.h5',
         'statepoint.07.h5', True)
    harness.main()
