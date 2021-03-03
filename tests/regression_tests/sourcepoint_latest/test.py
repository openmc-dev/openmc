import os
import sys

from tests.testing_harness import TestHarness


class SourcepointTestHarness(TestHarness):
    def _test_output_created(self):
        """Make sure statepoint.* and source* have been created."""
        TestHarness._test_output_created(self)
        source = glob.glob(os.path.join(os.getcwd(), 'source.*.h5'))
        assert len(source) == 1, 'Either multiple or no source files ' \
             'exist.'


def test_sourcepoint_latest():
    harness = TestHarness('statepoint.10.h5')
    harness.main()
