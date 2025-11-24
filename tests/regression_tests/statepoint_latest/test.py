import os
import glob

from tests.testing_harness import TestHarness


class StatepointLatestTestHarness(TestHarness):
    def _test_output_created(self):
        """Make sure running statepoint files have been created/rotated."""
        TestHarness._test_output_created(self)
        # Expect two rotating running files
        running = glob.glob(os.path.join(os.getcwd(), 'statepoint.running.*.h5'))
        assert len(running) == 2, 'Expected 2 running statepoint files to exist.'


def test_statepoint_latest():
    harness = TestHarness('statepoint.4.h5')
    harness.main()
