#!/usr/bin/env python

import glob
import os
import sys
sys.path.insert(0, os.pardir)
from testing_harness import TestHarness
from openmc.statepoint import StatePoint


class DistribDumpMatsTestHarness(TestHarness):
    def _test_output_created(self):
        """Make sure statepoint.* files have been created."""
        statepoint = glob.glob(os.path.join(os.getcwd(), self._sp_name))
        assert len(statepoint) == 1, 'Either multiple or no statepoint files ' \
             'exist.'
        assert statepoint[0].endswith('h5'), \
             'Statepoint file is not a HDF5 file.'
        if self._tallies:
            assert os.path.exists(os.path.join(os.getcwd(), 'tallies.out')), \
                 'Tally output file does not exist.'
		assert os.path.exists('material.m1.binary') or os.path.exists('material.m1.h5'), 'Materials file does not exist.'
		assert os.path.exists('material.m2.binary') or os.path.exists('material.m2.h5'), 'Materials file does not exist.'
		assert os.path.exists('material.m3.binary') or os.path.exists('material.m3.h5'), 'Materials file does not exist.'
	
	def _cleanup(self):
        """Delete statepoints, tally, and test files."""
        output = glob.glob(os.path.join(os.getcwd(), 'statepoint.*.*'))
		output += glob.glob(os.path.join(cwd, 'material.m*'))
        output.append(os.path.join(os.getcwd(), 'tallies.out'))
        output.append(os.path.join(os.getcwd(), 'results_test.dat'))
        for f in output:
            if os.path.exists(f):
                os.remove(f)


if __name__ == '__main__':
    harness = DistribDumpMatsTestHarness('statepoint.20.*')
    harness.main()
