#!/usr/bin/env python

import sys
sys.path.insert(0, '..')
from testing_harness import *


class StatepointRestartTestHarness(TestHarness):
    def execute_test(self):
        self._parse_args()
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

    def _run_openmc_restart(self):
        # Get the number of MPI processes.
        if self._opts.mpi_exec:
            mpi_procs = self._opts.mpi_np
        else:
            mpi_procs = 1

        # Get the name of the statepoint file.
        statepoint = glob.glob(os.path.join(os.getcwd(), self._sp_name))

        # Run OpenMC
        executor = Executor()
        returncode = executor.run_simulation(mpi_procs=mpi_procs,
                                             restart_file=statepoint)
        assert returncode == 0, 'OpenMC did not exit successfully.'


if __name__ == '__main__':
    harness = StatepointRestartTestHarness('statepoint.07.*', True)
    harness.execute_test()
