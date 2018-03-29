import glob
import os

from openmc import StatePoint

from tests.testing_harness import TestHarness


class EntropyTestHarness(TestHarness):
    def _get_results(self):
        """Digest info in the statepoint and return as a string."""
        # Read the statepoint file.
        statepoint = glob.glob(os.path.join(os.getcwd(), self._sp_name))[0]
        with StatePoint(statepoint) as sp:
            # Write out k-combined.
            outstr = 'k-combined:\n'
            outstr += '{:12.6E} {:12.6E}\n'.format(sp.k_combined.n, sp.k_combined.s)

            # Write out entropy data.
            outstr += 'entropy:\n'
            results = ['{:12.6E}'.format(x) for x in sp.entropy]
            outstr += '\n'.join(results) + '\n'

        return outstr


def test_entropy():
    harness = EntropyTestHarness('statepoint.10.h5')
    harness.main()
