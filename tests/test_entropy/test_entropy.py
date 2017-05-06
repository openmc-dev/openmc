#!/usr/bin/env python

import glob
import os
import sys
sys.path.insert(0, os.pardir)
from testing_harness import TestHarness
from openmc import StatePoint


class EntropyTestHarness(TestHarness):
    def _get_results(self):
        """Digest info in the statepoint and return as a string."""
        # Read the statepoint file.
        statepoint = glob.glob(os.path.join(os.getcwd(), self._sp_name))[0]
        with StatePoint(statepoint) as sp:
            # Write out k-combined.
            outstr = 'k-combined:\n'
            outstr += '{0:12.6E} {1:12.6E}\n'.format(*sp.k_combined)

            # Write out entropy data.
            outstr += 'entropy:\n'
            results = ['{0:12.6E}'.format(x) for x in sp.entropy]
            outstr += '\n'.join(results) + '\n'

        return outstr


if __name__ == '__main__':
    harness = EntropyTestHarness('statepoint.10.h5')
    harness.main()
