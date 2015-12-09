#!/usr/bin/env python

import glob
import os
import sys
sys.path.insert(0, os.pardir)
from testing_harness import TestHarness
import h5py

class DistribDumpMatsTestHarness(TestHarness):
    def _test_output_created(self):
        assert os.path.exists('materials-out.h5'), 'Materials file does not exist.'

    def _get_results(self):
        # set up output string
        outstr = ''
        # get material file contents
        matfile = 'materials-out.h5'
        file = h5py.File(matfile, 'r')
        for i in file.keys():
            group = file[i]
            outstr += '%s\n' % group.name
            outstr += '%s\n' % group['n_nuclides'].value
            outstr += '%s\n' % group['n_instances'].value
            for comp in group['comps']:
                outstr += ' '.join([str(comp)]) + '\n'        
        return outstr
                        
    def _cleanup(self):
        """Delete statepoints, tally, and test files."""
        output = glob.glob(os.path.join(os.getcwd(), 'statepoint.*.*'))
        output.append(os.path.join(os.getcwd(), 'materials-out.h5'))
        output.append(os.path.join(os.getcwd(), 'tallies.out'))
        output.append(os.path.join(os.getcwd(), 'results_test.dat'))
        for f in output:
            if os.path.exists(f):
                os.remove(f)


if __name__ == '__main__':
    harness = DistribDumpMatsTestHarness('statepoint.20.*')
    harness.main()
