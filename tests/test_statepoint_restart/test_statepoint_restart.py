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

            self._run_openmc_restart1()
            self._test_output_created()
            results = self._get_results()
            self._write_results(results)
            self._compare_results()

            self._run_openmc_restart2()
            self._test_output_created()
            results = self._get_results()
            self._write_results(results)
            self._compare_results()
        finally:
            self._cleanup()

    def _run_openmc_restart1(self):
        statepoint = glob.glob(os.path.join(os.getcwd(), self._sp_name))

        if self._opts.mpi_exec != '':
            proc = Popen([self._opts.mpi_exec, '-np', self._opts.mpi_np,
                          self._opts.exe, '-r', statepoint[0], os.getcwd()],
                         stderr=STDOUT, stdout=PIPE)
        else:
            proc = Popen([self._opts.exe, '-r', statepoint[0], os.getcwd()],
                         stderr=STDOUT, stdout=PIPE)
        print(proc.communicate()[0])
        returncode = proc.returncode
        assert returncode == 0, 'OpenMC did not exit successfully.'

    def _run_openmc_restart2(self):
        statepoint = glob.glob(os.path.join(os.getcwd(), self._sp_name))

        if self._opts.mpi_exec != '':
            proc = Popen([self._opts.mpi_exec, '-np', self._opts.mpi_np,
                          self._opts.exe, '--restart', statepoint[0],
                          os.getcwd()], stderr=STDOUT, stdout=PIPE)
        else:
            proc = Popen([self._opts.exe, '--restart', statepoint[0],
                          os.getcwd()], stderr=STDOUT, stdout=PIPE)
        print(proc.communicate()[0])
        returncode = proc.returncode
        assert returncode == 0, 'OpenMC did not exit successfully.'


if __name__ == '__main__':
    harness = StatepointRestartTestHarness('statepoint.07.*', True)
    harness.execute_test()
