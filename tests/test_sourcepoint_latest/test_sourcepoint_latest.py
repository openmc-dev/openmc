#!/usr/bin/env python

import sys
sys.path.insert(0, '..')
from testing_harness import *


class SourcepointTestHarness(TestHarness):
    def _test_output_created(self):
        """Make sure statepoint.* and source* have been created."""
        TestHarness._test_output_created(self)
        source = glob.glob(os.path.join(os.getcwd(), 'source.*'))
        assert len(source) == 1, 'Either multiple or no source files ' \
             'exist.'
        assert source[0].endswith('binary') \
             or source[0].endswith('h5'), \
             'Source file is not a binary or hdf5 file.'


if __name__ == '__main__':
    harness = TestHarness('statepoint.10.*')
    harness.main()
