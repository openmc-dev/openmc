#!/usr/bin/env python

import sys
sys.path.insert(0, '..')
from testing_harness import *


class OutputTestHarness(TestHarness):
    def _test_output_created(self):
        """Make sure output files have been created."""
        # Check for the statepoint.
        TestHarness._test_output_created(self)

        # Check for the summary.
        summary = glob.glob(os.path.join(os.getcwd(), 'summary.*'))
        assert len(summary) == 1, 'Either multiple or no summary file exists.'
        assert summary[0].endswith('out') or summary[0].endswith('h5'),\
            'Summary file is not a binary or hdf5 file.'

        # Check for the cross sections.
        assert os.path.exists(os.path.join(os.getcwd(), 'cross_sections.out')),\
            'Cross section output file does not exist.'

    def _cleanup(self):
        TestHarness._cleanup(self)
        output = glob.glob(os.path.join(os.getcwd(), 'summary.*'))
        output.append(os.path.join(os.getcwd(), 'cross_sections.out'))
        for f in output:
            if os.path.exists(f):
                os.remove(f)


if __name__ == '__main__':
    harness = OutputTestHarness('statepoint.10.*')
    harness.main()
