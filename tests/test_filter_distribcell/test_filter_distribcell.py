#!/usr/bin/env python

import hashlib
import sys

sys.path.insert(0, '..')
from testing_harness import *


class DistribcellTestHarness(TestHarness):
    def __init__(self):
        self._sp_name = None
        self._tallies = True
        self._opts = None
        self._args = None


    def execute_test(self):
        self._parse_args()
        base_dir = os.getcwd()
        try:
            self._run_openmc()
            os.chdir(base_dir)
            self._test_output_created()
            os.chdir(base_dir)
            self._get_results()
            os.chdir(base_dir)
            self._compare_results()
        finally:
            os.chdir(base_dir)
            self._cleanup()


    def _run_openmc(self):
        dirs = ('case-1', '../case-2', '../case-3')
        for d in dirs:
            os.chdir(d)
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
        """Make sure statepoint files have been created."""
        dirs = ('case-1', '../case-2', '../case-3')
        sps = ('statepoint.1.*', 'statepoint.1.*', 'statepoint.3.*')
        tallies_present = (True, True, False)
        for i in range(len(dirs)):
            os.chdir(dirs[i])
            self._tallies = tallies_present[i]
            self._sp_name = sps[i]
            TestHarness._test_output_created(self)
        self._tallies = True


    def _get_results(self):
        dirs = ('case-1', '../case-2', '../case-3')
        sps = ('statepoint.1.*', 'statepoint.1.*', 'statepoint.3.*')
        for i in range(len(dirs)):
            os.chdir(dirs[i])
            self._sp_name = sps[i]

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

            if i == 2:
                sha512 = hashlib.sha512()
                sha512.update(outstr.encode('utf-8'))
                outstr = sha512.hexdigest()

            # Write results to a file.
            with open('results_test.dat','w') as fh:
                fh.write(outstr)


    def _compare_results(self):
        dirs = ('case-1', '../case-2', '../case-3')
        for d in dirs:
            os.chdir(d)
            TestHarness._compare_results(self)

    def _cleanup(self):
        dirs = ('case-1', '../case-2', '../case-3')
        for d in dirs:
            os.chdir(d)
            TestHarness._cleanup(self)


if __name__ == '__main__':
    harness = DistribcellTestHarness()
    harness.execute_test()
