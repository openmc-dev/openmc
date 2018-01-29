#!/usr/bin/env python

import glob
import os

from tests.testing_harness import *


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
        """Run OpenMC with the appropriate arguments and check the outputs."""
        try:
            self._run_openmc()
            self._test_output_created()
            self._run_openmc_restart()
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
            self._run_openmc_restart()
            results = self._get_results()
            self._write_results(results)
            self._overwrite_results()
        finally:
            self._cleanup()

    def _test_output_created(self):
        """Make sure statepoint and source files have been created."""
        statepoint = glob.glob(os.path.join(os.getcwd(), self._sp_name))
        assert len(statepoint) == 1, 'Either multiple or no statepoint files ' \
             'exist.'
        assert statepoint[0].endswith('h5'), \
             'Statepoint file is not a HDF5 file.'

        source = glob.glob(os.path.join(os.getcwd(), 'source.10.*'))
        assert len(source) == 1, 'Either multiple or no source files exist.'
        assert source[0].endswith('h5'), \
             'Source file is not a HDF5 file.'

    def _run_openmc_restart(self):
        # Get the name of the source file.
        source = glob.glob(os.path.join(os.getcwd(), 'source.10.*'))

        # Write the new settings.xml file.
        with open('settings.xml','w') as fh:
            fh.write(settings2.format(source[0].split('.')[-1]))

        # Run OpenMC.
        self._run_openmc()

    def _cleanup(self):
        TestHarness._cleanup(self)
        output = glob.glob(os.path.join(os.getcwd(), 'source.*'))
        for f in output:
            if os.path.exists(f):
                os.remove(f)
        with open('settings.xml','w') as fh:
            fh.write(settings1)


def test_source_file():
    harness = SourceFileTestHarness('statepoint.10.h5')
    harness.main()
