#!/usr/bin/env python

import os
import sys
import glob
import shutil

sys.path.insert(0, os.pardir)
from testing_harness import TestHarness
import openmc
import numpy as np


class DomainDecomTestHarness(TestHarness):
    def execute_test(self):
        """Run OpenMC with the appropriate arguments and check the outputs."""
        if self._opts.mpi_exec is None:
            # This test is only relevant if MPI is enabled
            return
        base_dir = os.getcwd()
        try:
            # Run the 1-domain case
            os.chdir('1_domain')
            shutil.copyfile("../results_true.dat", "./results_true.dat")
            self._run_openmc(1)
            self._test_output_created(1)
            results = self._get_results()
            self._write_results(results)
            self._compare_results()

            # Run the 4-domain case
            os.chdir('../4_domains')
            shutil.copyfile("../results_true.dat", "./results_true.dat")
            self._run_openmc(4)
            self._test_output_created(4)
            results = self._get_results()
            self._write_results(results)
            self._compare_results()

        finally:
            os.chdir(base_dir)
            os.chdir('1_domain')
            self._cleanup()
            os.chdir('../4_domains')
            self._cleanup()

    def update_results(self):
        """Update the results_true using the current version of OpenMC."""
        if self._opts.mpi_exec is None:
            # This test is only relevant if MPI is enabled
            return
        base_dir = os.getcwd()
        try:
            # Run the 1-domain case
            os.chdir('1_domain')
            self._run_openmc(1)
            self._test_output_created(1)
            results = self._get_results()
            self._write_results(results)
            self._overwrite_results()
            shutil.copyfile("./results_true.dat", "../results_true.dat")

        finally:
            os.chdir(base_dir)
            os.chdir('1_domain')
            self._cleanup()

    def _run_openmc(self, domain_num):
        if domain_num == 1:
            p_num = 3
        elif domain_num == 4:
            p_num = 5

        returncode = openmc.run(mpi_procs=p_num, threads=1,
                                openmc_exec=self._opts.exe,
                                mpi_exec=self._opts.mpi_exec)
        assert returncode == 0, 'OpenMC did not exit successfully.'

    def _test_output_created(self, domain_num):
        """Make sure statepoint.* and tallies.out have been created."""
        statepoint = glob.glob(os.path.join(os.getcwd(), self._sp_name))
        tallies = glob.glob(os.path.join(os.getcwd(), 'tallies.domain_*'))
        assert len(statepoint) == 1, \
            'Either multiple or no statepoint files exist'
        assert len(tallies) == domain_num, \
            'Wrong number of tally out files: %s' % len(tallies)
        assert statepoint[0].endswith('h5'), \
            'Statepoint file is not a HDF5 file.'

    def _get_results(self):
        """Digest info in the statepoint and return as a string."""
        # Read the statepoint file
        statepoint = glob.glob(os.path.join(os.getcwd(), self._sp_name))[0]
        sp = openmc.StatePoint(statepoint)

        # Write out k-combined
        outstr = 'k-combined:\n'
        form = '{0:12.6E} {1:12.6E}\n'
        outstr += form.format(sp.k_combined[0], sp.k_combined[1])

        # Write out tally data
        # Std_dev will not be written as it is not reproducible
        if self._tallies:
            tally_num = 1
            for tally_ind in sp.tallies:
                tally = sp.tallies[tally_ind]
                results = tally.sum.ravel()
                results = ['{0:12.6E}'.format(x) for x in results]

                outstr += 'tally ' + str(tally_num) + ':\n'
                outstr += '\n'.join(results) + '\n'
                tally_num += 1

        return outstr

    def _cleanup(self):
        """Delete statepoints, tally, and test files."""
        output = glob.glob(os.path.join(os.getcwd(), 'statepoint.*.*'))
        output.append(os.path.join(os.getcwd(), 'summary.h5'))
        output.append(os.path.join(os.getcwd(), 'results*.dat'))
        output += glob.glob(os.path.join(os.getcwd(), 'tallies*.out'))
        output += glob.glob(os.path.join(os.getcwd(), 'volume_*.h5'))
        for f in output:
            if os.path.exists(f):
                os.remove(f)

if __name__ == '__main__':
    harness = DomainDecomTestHarness('statepoint.20.*', True)
    harness.main()