#!/usr/bin/env python

import sys
sys.path.insert(0, '..')
from testing_harness import *


settings1="""<?xml version="1.0"?>
<settings>
  <state_point batches="10" />
  <source_point separate="true" />
  <eigenvalue>
    <batches>10</batches>
    <inactive>5</inactive>
    <particles>1000</particles>
  </eigenvalue>
  <source>
    <space type="box">
      <parameters>-4 -4 -4  4  4  4</parameters>
    </space>
  </source>
</settings>
"""

settings2 = """<?xml version="1.0"?>
<settings>
  <eigenvalue>
    <batches>10</batches>
    <inactive>5</inactive>
    <particles>1000</particles>
  </eigenvalue>
  <source>
    <file> source.10.{0} </file>
  </source>
</settings>
"""


class SourceFileTestHarness(TestHarness):
    def execute_test(self):
        self._parse_args()
        try:
            self._run_openmc()
            self._test_output_created()
            self._run_openmc_restart()
            results = self._get_results()
            self._write_results(results)
            self._compare_results()
        finally:
            self._cleanup()

    def _test_output_created(self):
        """Make sure statepoint and source files have been created."""
        statepoint = glob.glob(os.path.join(os.getcwd(), self._sp_name))
        assert len(statepoint) == 1, 'Either multiple or no statepoint files ' \
             'exist.'
        assert statepoint[0].endswith('binary') \
             or statepoint[0].endswith('h5'), \
             'Statepoint file is not a binary or hdf5 file.'

        source = glob.glob(os.path.join(os.getcwd(), 'source.10.*'))
        assert len(source) == 1, 'Either multiple or no source files exist.'
        assert source[0].endswith('binary') \
             or source[0].endswith('h5'), \
             'Source file is not a binary or hdf5 file.'

    def _run_openmc_restart(self):
        # Get the number of MPI processes.
        if self._opts.mpi_exec:
            mpi_procs = self._opts.mpi_np
        else:
            mpi_procs = 1

        # Get the name of the source file.
        source = glob.glob(os.path.join(os.getcwd(), 'source.10.*'))

        # Write the new settings.xml file.
        with open('settings.xml','w') as fh:
            fh.write(settings2.format(source[0].split('.')[-1]))

        # Run OpenMC.
        executor = Executor()
        returncode = executor.run_simulation(mpi_procs=mpi_procs)
        assert returncode == 0, 'OpenMC did not exit successfully.'

    def _cleanup(self):
        TestHarness._cleanup(self)
        output = glob.glob(os.path.join(os.getcwd(), 'source.*'))
        for f in output:
            if os.path.exists(f):
                os.remove(f)
        with open('settings.xml','w') as fh:
            fh.write(settings1)


if __name__ == '__main__':
    harness = SourceFileTestHarness('statepoint.10.*')
    harness.execute_test()
