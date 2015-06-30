#!/usr/bin/env python 
import sys
sys.path.insert(0, '..')
from testing_harness import *


class StatepointTestHarness(TestHarness):
    def __init__(self):
        self._sp_name = None
        self._tallies = False
        self._opts = None
        self._args = None

    def _test_output_created(self):
        """Make sure statepoint files have been created."""
        sps = ('statepoint.03.*', 'statepoint.06.*', 'statepoint.09.*')
        for sp in sps:
            self._sp_name = sp
            TestHarness._test_output_created(self)


if __name__ == '__main__':
    harness = StatepointTestHarness()
    harness.main()
