import glob

from openmc import StatePoint

from tests.testing_harness import TestHarness


class SourcepointTestHarness(TestHarness):
    def _test_output_created(self):
        """Make sure statepoint files have been created."""
        statepoint = glob.glob('statepoint.*.h5')
        assert len(statepoint) == 5, 'Five statepoint files must exist.'

    def _get_results(self):
        """Digest info in the statepoint and return as a string."""
        # Get the eigenvalue information.
        outstr = TestHarness._get_results(self)

        # Read the statepoint file.
        with StatePoint(self._sp_name) as sp:
            # Add the source information.
            xyz = sp.source[0]['r']
            outstr += ' '.join(['{0:12.6E}'.format(x) for x in xyz])
            outstr += "\n"

        return outstr


def test_sourcepoint_batch():
    harness = SourcepointTestHarness('statepoint.08.h5')
    harness.main()
