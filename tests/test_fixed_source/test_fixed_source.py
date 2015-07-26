#!/usr/bin/env python

import sys
sys.path.insert(0, '..')
from testing_harness import *


class FixedSourceTestHarness(TestHarness):
    def _get_results(self):
        """Digest info in the statepoint and return as a string."""
        # Read the statepoint file.
        statepoint = glob.glob(os.path.join(os.getcwd(), self._sp_name))[0]
        sp = StatePoint(statepoint)
        sp.read_results()

        # Write out tally data.
        outstr = ''
        if self._tallies:
            tally_num = 1
            for tally_ind in sp._tallies:
                tally = sp._tallies[tally_ind]
                results = np.zeros((tally._sum.size*2, ))
                results[0::2] = tally._sum.ravel()
                results[1::2] = tally._sum_sq.ravel()
                results = ['{0:12.6E}'.format(x) for x in results]

                outstr += 'tally ' + str(tally_num) + ':\n'
                outstr += '\n'.join(results) + '\n'
                tally_num += 1

        return outstr


if __name__ == '__main__':
    harness = FixedSourceTestHarness('statepoint.10.*', True)
    harness.main()
