#!/usr/bin/env python

import glob
import os
import sys

sys.path.insert(0, '..')
from testing_harness import *


class FixedSourceTestHarness(TestHarness):
    def _get_results(self):
        """Digest info in the statepoint and return as a string."""
        # Read the statepoint file.
        statepoint = glob.glob(os.path.join(os.getcwd(), self._sp_name))[0]
        sp = StatePoint(statepoint)

        # Write out tally data.
        outstr = ''
        if self._tallies:
            tally_num = 1
            for tally_ind in sp.tallies:
                tally = sp.tallies[tally_ind]
                results = np.zeros((tally.sum.size*2, ))
                results[0::2] = tally.sum.ravel()
                results[1::2] = tally.sum_sq.ravel()
                results = ['{0:12.6E}'.format(x) for x in results]

                outstr += 'tally ' + str(tally_num) + ':\n'
                outstr += '\n'.join(results) + '\n'
                tally_num += 1

        outstr += 'leakage:\n'
        outstr += '{0:12.6E}'.format(sp.global_tallies[3][0]) + '\n'
        outstr += '{0:12.6E}'.format(sp.global_tallies[3][1]) + '\n'

        return outstr


if __name__ == '__main__':
    harness = FixedSourceTestHarness('statepoint.10.*', True)
    harness.main()
