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
        self._parse_args()
        base_dir = os.getcwd()
        try:
            dirs = ('case-1', '../case-2', '../case-3')
            sps = ('statepoint.1.*', 'statepoint.1.*', 'statepoint.3.*')
            tallies_present = (True, True, False)
            hash_out = (False, False, True)
            for i in range(len(dirs)):
                os.chdir(dirs[i])
                self._sp_name = sps[i]
                self._tallies = tallies_present[i]

                self._run_openmc()
                self._test_output_created()
                results = self._get_results(hash_out[i])
                self._write_results(results)
                self._compare_results()
        finally:
            os.chdir(base_dir)
            for i in range(len(dirs)):
                os.chdir(dirs[i])
                self._cleanup()


if __name__ == '__main__':
    harness = DistribcellTestHarness()
    harness.execute_test()
