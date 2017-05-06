#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
from testing_harness import TestHarness


class SourcepointTestHarness(TestHarness):
    def _test_output_created(self):
        """Make sure statepoint.* and source* have been created."""
        TestHarness._test_output_created(self)
        source = glob.glob(os.path.join(os.getcwd(), 'source.*'))
        assert len(source) == 1, 'Either multiple or no source files ' \
             'exist.'
        assert source[0].endswith('h5'), \
             'Source file is not a HDF5 file.'


if __name__ == '__main__':
    harness = TestHarness('statepoint.10.h5')
    harness.main()
