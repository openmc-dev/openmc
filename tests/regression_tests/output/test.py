import os
import glob

from tests.testing_harness import TestHarness


class OutputTestHarness(TestHarness):
    def _test_output_created(self):
        """Make sure output files have been created."""
        # Check for the statepoint.
        TestHarness._test_output_created(self)

        # Check for the summary.
        summary = glob.glob(os.path.join(os.getcwd(), 'summary.*'))
        assert len(summary) == 1, 'Either multiple or no summary file exists.'
        assert summary[0].endswith('h5'),\
            'Summary file is not a HDF5 file.'

    def _cleanup(self):
        TestHarness._cleanup(self)
        f = 'summary.h5'
        if os.path.exists(f):
            os.remove(f)


def test_output():
    harness = OutputTestHarness('statepoint.10.h5')
    harness.main()
