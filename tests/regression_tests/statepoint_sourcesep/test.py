import glob
import os

from tests.testing_harness import TestHarness


class SourcepointTestHarness(TestHarness):
    def _test_output_created(self):
        """Make sure statepoint.* and source* have been created."""
        TestHarness._test_output_created(self)
        source = glob.glob('source.*.h5')
        assert len(source) == 1, 'Either multiple or no source files ' \
             'exist.'

    def _cleanup(self):
        TestHarness._cleanup(self)
        output = glob.glob('source.*.h5')
        for f in output:
            if os.path.exists(f):
                os.remove(f)


def test_statepoint_sourcesep():
    harness = SourcepointTestHarness('statepoint.10.h5')
    harness.main()
