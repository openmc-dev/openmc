#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
from testing_harness import TestHarness


class StatepointTestHarness(TestHarness):
    def __init__(self):
        super(StatepointTestHarness, self).__init__(None)

    def _test_output_created(self):
        """Make sure statepoint files have been created."""
        sps = ('statepoint.03.h5', 'statepoint.06.h5', 'statepoint.09.h5')
        for sp in sps:
            self._sp_name = sp
            TestHarness._test_output_created(self)


if __name__ == '__main__':
    harness = StatepointTestHarness()
    harness.main()
