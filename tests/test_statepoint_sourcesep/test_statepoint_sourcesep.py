#!/usr/bin/env python

import glob
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

    def _cleanup(self):
        TestHarness._cleanup(self)
        output = glob.glob(os.path.join(os.getcwd(), 'source.*'))
        for f in output:
            if os.path.exists(f):
                os.remove(f)


if __name__ == '__main__':
    harness = SourcepointTestHarness('statepoint.10.h5')
    harness.main()
