#!/usr/bin/env python

import sys
sys.path.insert(0, '..')
from testing_harness import *


class CMFDTestHarness(TestHarness):
    def _get_results(self):
        """Digest info in the statepoint and create a simpler ASCII file."""
        # Read the statepoint file.
        statepoint = glob.glob(os.path.join(os.getcwd(), self._sp_name))[0]
        sp = StatePoint(statepoint)
        sp.read_results()

        # Write out k-combined.
        outstr = 'k-combined:\n'
        form = '{0:12.6E} {1:12.6E}\n'
        outstr += form.format(sp.k_combined[0], sp.k_combined[1])

        # Write out tally data.
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

        # Write out CMFD data.
        outstr += 'cmfd indices\n'
        outstr += '\n'.join(['{0:12.6E}'.format(x) for x in sp._cmfd_indices])
        outstr += '\nk cmfd\n'
        outstr += '\n'.join(['{0:12.6E}'.format(x) for x in sp._k_cmfd])
        outstr += '\ncmfd entropy\n'
        outstr += '\n'.join(['{0:12.6E}'.format(x) for x in sp._cmfd_entropy])
        outstr += '\ncmfd balance\n'
        outstr += '\n'.join(['{0:12.6E}'.format(x) for x in sp._cmfd_balance])
        outstr += '\ncmfd dominance ratio\n'
        outstr += '\n'.join(['{0:10.3E}'.format(x) for x in sp._cmfd_dominance])
        outstr += '\ncmfd openmc source comparison\n'
        outstr += '\n'.join(['{0:12.6E}'.format(x) for x in sp._cmfd_srccmp])
        outstr += '\ncmfd source\n'
        cmfdsrc = np.reshape(sp._cmfd_src, np.product(sp._cmfd_indices),
                             order='F')
        outstr += '\n'.join(['{0:12.6E}'.format(x) for x in cmfdsrc])
        outstr += '\n'

        # Write results to a file.
        with open('results_test.dat','w') as fh:
            fh.write(outstr)


if __name__ == '__main__':
    harness = CMFDTestHarness('statepoint.20.*', True)
    harness.execute_test()
