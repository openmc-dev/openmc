#!/usr/bin/env python

import os
import sys
import glob

sys.path.insert(0, os.pardir)
from testing_harness import *
from openmc.executor import Executor
from openmc.statepoint import StatePoint
import numpy as np

def order_by(arr, ordering):
    ordered = np.zeros(arr.shape)
    for i,val in enumerate(arr):
        ordered[ordering[i]-1] = val
    return ordered 

class DomainDecomOTFTalliesTestHarness(TestHarness):
    def __init__(self, statepoint_name):
        super(DomainDecomOTFTalliesTestHarness, self).__init__(statepoint_name, True)
        
    def execute_test(self):
        """Run OpenMC with the appropriate arguments and check the outputs."""
        if self._opts.mpi_exec is None:
            # This test is only relevant if MPI is enabled
            return
        base_dir = os.getcwd()
        try:
            # Run the 1-domain case
            os.chdir('1_domain')
            self._run_openmc(1)
            self._test_output_created(1)
            results = self._get_results(1)
            self._write_results(results)
            self._compare_results()
            
            # Run the 4-domain case
            os.chdir('../4_domains')
            self._run_openmc(4)
            self._test_output_created(4)
            results = self._get_results(4)
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
            results = self._get_results(1)
            self._write_results(results)
            self._overwrite_results()
            
            # Run the 4-domain case
            os.chdir('../4_domains')
            self._run_openmc(4)
            self._test_output_created(4)
            results = self._get_results(4)
            self._write_results(results)
            self._overwrite_results()
            
        finally:
            os.chdir(base_dir)
            os.chdir('1_domain')
            self._cleanup()
            os.chdir('../4_domains')
            self._cleanup()

    def _run_openmc(self, domain_num):
        if domain_num == 1:
            p_num = 3
        elif domain_num == 4:
            p_num = 5
        
        executor = Executor()
        returncode = executor.run_simulation(mpi_procs=p_num,
                                             openmc_exec=self._opts.exe,
                                             mpi_exec=self._opts.mpi_exec)
        assert returncode == 0, 'OpenMC did not exit successfully.'

    def _test_output_created(self, domain_num):
        """Make sure statepoint.* and tallies.out have been created."""
        statepoint = glob.glob(os.path.join(os.getcwd(), self._sp_name))
        if domain_num == 1:
            assert len(statepoint) == 1, 'Either multiple or no statepoint files exist.'
            assert statepoint[0].endswith('h5'), 'Statepoint file is not a HDF5 file.'
            assert os.path.exists(os.path.join(os.getcwd(), 'tallies.domain_1.out')),'Tally output file does not exist.'
        elif domain_num == 4:
            assert len(statepoint) == 4, 'Wrong number of statepoint files: %s' % \
                len(statepoint)
            assert statepoint[0].endswith('h5'), 'Statepoint file is not a HDF5 file.'
            assert os.path.exists(os.path.join(os.getcwd(), 'tallies.domain_1.out')), 'Tally output file does not exist.'
            assert os.path.exists(os.path.join(os.getcwd(), 'tallies.domain_2.out')), 'Tally output file does not exist.'
            assert os.path.exists(os.path.join(os.getcwd(), 'tallies.domain_3.out')), 'Tally output file does not exist.'
            assert os.path.exists(os.path.join(os.getcwd(), 'tallies.domain_4.out')), 'Tally output file does not exist.'
    
    def _get_results(self, domain_num):
        """Digest info in the statepoint and return as a string."""
        if domain_num == 1:
            # Read the statepoint file.
            statepoint = glob.glob(os.path.join(os.getcwd(), 'statepoint.20.domain_1.h5'))[0]
            sp = StatePoint(statepoint)
            # extract tally results (means only) and convert to vector
            results = sp.tallies[1].mean[:,:,0]
            otf_filter_bin_map = sp._f['tallies/tally 1/otf_filter_bin_map']
            results = order_by(results, otf_filter_bin_map)
            shape = results.shape
            size = (np.product(shape))
            results = np.reshape(results, size)
            # set up output string
            outstr = ''
            # write out k-combined
            outstr += 'k-combined:\n'
            outstr += "{0:12.6E} {1:12.6E}\n".format(sp.k_combined[0], sp.k_combined[1])
            # write out tally results
            outstr += 'tallies:\n'
            for item in results:
              outstr += "{0:12.6E}\n".format(item)
            return outstr
        elif domain_num == 4:
            # read in statepoint files
            spfile = 'statepoint.20.domain_1.h5'
            statepoint = glob.glob(os.path.join(os.getcwd(), spfile))[0]
            sp1 = StatePoint(statepoint)
            spfile = 'statepoint.20.domain_2.h5'
            statepoint = glob.glob(os.path.join(os.getcwd(), spfile))[0]
            sp2 = StatePoint(statepoint)
            spfile = 'statepoint.20.domain_3.h5'
            statepoint = glob.glob(os.path.join(os.getcwd(), spfile))[0]
            sp3 = StatePoint(statepoint)
            spfile = 'statepoint.20.domain_4.h5'
            statepoint = glob.glob(os.path.join(os.getcwd(), spfile))[0]
            sp4 = StatePoint(statepoint)
            # extract tally results, means only since sum_sq won't match the 1_domain case
            mean1=sp1._f['tallies/tally 1/results'].value['sum']/sp1.tallies[1].num_realizations
            mean2=sp2._f['tallies/tally 1/results'].value['sum']/sp2.tallies[1].num_realizations
            mean3=sp3._f['tallies/tally 1/results'].value['sum']/sp3.tallies[1].num_realizations
            mean4=sp4._f['tallies/tally 1/results'].value['sum']/sp4.tallies[1].num_realizations
            results = [mean1, mean2, mean3, mean4]
            # combine results for bins that were on more than one domain, in real_bin order
            otf_filter_bin_map1 = sp1._f['tallies/tally 1/otf_filter_bin_map']
            otf_filter_bin_map2 = sp2._f['tallies/tally 1/otf_filter_bin_map']
            otf_filter_bin_map3 = sp3._f['tallies/tally 1/otf_filter_bin_map']
            otf_filter_bin_map4 = sp4._f['tallies/tally 1/otf_filter_bin_map']
            maps = np.array([otf_filter_bin_map1,
                             otf_filter_bin_map2,
                             otf_filter_bin_map3,
                             otf_filter_bin_map4])
            maxbin = max(maps.flatten())
            # maxbin = max(max(maps))
            tmp_results = np.zeros((maxbin, 1))
            for i in range(maxbin):
                for j, map_ in enumerate(maps):
                    for k, index in enumerate(map_):
                        if i+1 == index:
                            tmp_results[i] += results[j][k][0]
            results = np.reshape(tmp_results, (np.product(tmp_results.shape)))
            # set up output string
            outstr = ''
            # write out k-combined
            outstr += 'k-combined:\n'
            outstr += "{0:12.6E} {1:12.6E}\n".format(sp1.k_combined[0], sp1.k_combined[1])
            # write out tally results
            outstr += 'tallies:\n'
            for item in results:
                outstr += "{0:12.6E}\n".format(item)
            return outstr

    def _cleanup(self):
        """Delete statepoints, tally, and test files."""
        output = glob.glob(os.path.join(os.getcwd(), 'statepoint.*.*'))
        output.append(os.path.join(os.getcwd(), 'tallies.domain*'))
        output.append(os.path.join(os.getcwd(), 'results_test.dat'))
        for f in output:
            if os.path.exists(f):
                os.remove(f)

if __name__ == '__main__':
    harness = DomainDecomOTFTalliesTestHarness('statepoint.20.*')
    harness.main()