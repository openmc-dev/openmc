#!/usr/bin/env python

import sys
sys.path.insert(0, '..')
from testing_harness import *


class PlotTestHarness(TestHarness):
    def execute_test(self):
        self._parse_args()
        try:
            self._run_openmc()
            self._test_output_created()
        finally:
            self._cleanup()


    def _run_openmc(self):
        if self._opts.mpi_exec != '':
            proc = Popen([self._opts.mpi_exec, '-np', self._opts.mpi_np,
                          self._opts.exe, '-p', os.getcwd()],
                         stderr=STDOUT, stdout=PIPE)
        else:
            proc = Popen([self._opts.exe, '-p', os.getcwd()],
                         stderr=STDOUT, stdout=PIPE)
        print(proc.communicate()[0])
        returncode = proc.returncode
        assert returncode == 0, 'OpenMC did not exit successfully.'


    def _test_output_created(self):
        """Make sure *.ppm has been created."""
        assert os.path.exists(os.path.join(os.getcwd(), '1_plot.ppm')), \
             'Plot output file does not exist.'


    def _cleanup(self):
        TestHarness._cleanup(self)
        output = glob.glob(os.path.join(os.getcwd(), '*.ppm'))
        for f in output:
            if os.path.exists(f):
                os.remove(f)


if __name__ == '__main__':
    harness = PlotTestHarness('statepoint.10.*', True)
    harness.execute_test()
