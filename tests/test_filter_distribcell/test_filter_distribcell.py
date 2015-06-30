#!/usr/bin/env python

import hashlib
import sys

sys.path.insert(0, '..')
from testing_harness import *


class DistribcellTestHarness(TestHarness):
    def __init__(self):
        self._sp_name = None
        self._tallies = True
        self._opts = None
        self._args = None

    def execute_test(self):
        """Run OpenMC with the appropriate arguments and check the outputs."""
        base_dir = os.getcwd()
        try:
            dirs = ('case-1', '../case-2', '../case-3')
            sps = ('statepoint.1.*', 'statepoint.1.*', 'statepoint.3.*')
            tallies_out_present = (True, True, False)
            hash_out = (False, False, True)
            for i in range(len(dirs)):
                os.chdir(dirs[i])
                self._sp_name = sps[i]

                self._run_openmc()
                self._test_output_created(tallies_out_present[i])
                results = self._get_results(hash_out[i])
                self._write_results(results)
                self._compare_results()
        finally:
            os.chdir(base_dir)
            for i in range(len(dirs)):
                os.chdir(dirs[i])
                self._cleanup()

    def update_results(self):
        """Update the results_true using the current version of OpenMC."""
        base_dir = os.getcwd()
        try:
            dirs = ('case-1', '../case-2', '../case-3')
            sps = ('statepoint.1.*', 'statepoint.1.*', 'statepoint.3.*')
            tallies_out_present = (True, True, False)
            hash_out = (False, False, True)
            for i in range(len(dirs)):
                os.chdir(dirs[i])
                self._sp_name = sps[i]

                self._run_openmc()
                self._test_output_created(tallies_out_present[i])
                results = self._get_results(hash_out[i])
                self._write_results(results)
                self._overwrite_results()
        finally:
            os.chdir(base_dir)
            for i in range(len(dirs)):
                os.chdir(dirs[i])
                self._cleanup()

    def _test_output_created(self, tallies_out_present):
        """Make sure statepoint.* and tallies.out have been created."""
        statepoint = glob.glob(os.path.join(os.getcwd(), self._sp_name))
        assert len(statepoint) == 1, 'Either multiple or no statepoint files ' \
             'exist.'
        assert statepoint[0].endswith('binary') \
             or statepoint[0].endswith('h5'), \
             'Statepoint file is not a binary or hdf5 file.'
        if tallies_out_present:
            assert os.path.exists(os.path.join(os.getcwd(), 'tallies.out')), \
                 'Tally output file does not exist.'


if __name__ == '__main__':
    harness = DistribcellTestHarness()
    harness.main()
