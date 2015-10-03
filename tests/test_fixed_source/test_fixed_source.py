#!/usr/bin/env python

import glob
import os
import sys
import numpy as np
sys.path.insert(0, os.pardir)
from testing_harness import TestHarness
from openmc.statepoint import StatePoint


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

        gt = sp.global_tallies
        outstr += 'leakage:\n'
        outstr += '{0:12.6E}'.format(gt[gt['name'] == b'leakage'][0]['sum']) + '\n'
        outstr += '{0:12.6E}'.format(gt[gt['name'] == b'leakage'][0]['sum_sq']) + '\n'

        return outstr


if __name__ == '__main__':
    harness = FixedSourceTestHarness('statepoint.10.*', True)
    harness.main()
