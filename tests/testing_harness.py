from __future__ import print_function

import filecmp
import glob
from optparse import OptionParser
import os
from subprocess import Popen, STDOUT, PIPE, call
import sys

import numpy as np

sys.path.insert(0, '../..')
from openmc.statepoint import StatePoint


class TestHarness():
    def __init__(self, statepoint_name, tallies_present=False):
        self._sp_name = statepoint_name
        self._tallies = tallies_present
        self._opts = None
        self._args = None


    def execute_test(self):
        self._parse_args()
        try:
            self._run_openmc()
            self._test_output_created()
            self._get_results()
            self._compare_results()
        finally:
            self._cleanup()


    def _parse_args(self):
        parser = OptionParser()
        parser.add_option('--mpi_exec', dest='mpi_exec', default='')
        parser.add_option('--mpi_np', dest='mpi_np', default='3')
        parser.add_option('--exe', dest='exe')
        (self._opts, self._args) = parser.parse_args()
        if self._opts.exe is None:
            raise Exception('Must specify OpenMC executable from command line '
                            'with --exe.')


    def _run_openmc(self):
        if self._opts.mpi_exec != '':
            proc = Popen([self._opts.mpi_exec, '-np', self._opts.mpi_np,
                          self._opts.exe, os.getcwd()],
                         stderr=STDOUT, stdout=PIPE)
        else:
            proc = Popen([self._opts.exe, os.getcwd()],
                         stderr=STDOUT, stdout=PIPE)
        print(proc.communicate()[0])
        returncode = proc.returncode
        assert returncode == 0, 'OpenMC did not exit successfully.'


    def _test_output_created(self):
        """Make sure statepoint.* and tallies.out have been created."""
        statepoint = glob.glob(os.path.join(os.getcwd(), self._sp_name))
        assert len(statepoint) == 1, 'Either multiple or no statepoint files ' \
             'exist.'
        assert statepoint[0].endswith('binary') \
             or statepoint[0].endswith('h5'), \
             'Statepoint file is not a binary or hdf5 file.'
        if self._tallies:
            assert os.path.exists(os.path.join(os.getcwd(), 'tallies.out')), \
                 'Tally output file does not exist.'


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

        # Write results to a file.
        with open('results_test.dat','w') as fh:
            fh.write(outstr)


    def _compare_results(self):
        compare = filecmp.cmp('results_test.dat', 'results_true.dat')
        if not compare:
            os.rename('results_test.dat', 'results_error.dat')
        assert compare, 'Results do not agree.'


    def _cleanup(self):
        output = glob.glob(os.path.join(os.getcwd(), 'statepoint.*.*'))
        output.append(os.path.join(os.getcwd(), 'tallies.out'))
        output.append(os.path.join(os.getcwd(), 'results_test.dat'))
        for f in output:
            if os.path.exists(f):
                os.remove(f)
