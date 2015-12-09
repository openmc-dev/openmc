#!/usr/bin/env python

import os
import sys
import glob

sys.path.insert(0, os.pardir)
from testing_harness import *
from openmc.executor import Executor
from openmc.statepoint import StatePoint
import numpy as np


class DomainDecomOTFTalliesTestHarness(TestHarness):
    def __init__(self, statepoint_name):
        self._sp_name = statepoint_name
        self._tallies = True
        self._opts = None
        self._args = None
    
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
            results = self._get_results()
            self._write_results(results)
            self._compare_results()
            
            # Run the 4-domain case
            os.chdir('../4_domains')
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
            
            # Run the 4-domain case
            os.chdir('../4_domains')
            self._run_openmc(4)
            self._test_output_created(4)
            results = self._get_results()
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
    
    def _get_results(self):
        """Digest info in the statepoint and return as a string."""
        # Read the statepoint file.
        statepoint = glob.glob(os.path.join(os.getcwd(), 'statepoint.20.domain_1.h5'))[0]
        sp = StatePoint(statepoint)
        sp.read_domain_tally_metadata()
        #sp.read_results()
        # set up output string
        outstr = ''
        spec_list = [('distribcell', [0, 1, (6,1,1,1), 5, 2, (4,1,2,1), 1, 101])]
        outstr += '%s\n' % sp.get_value(1, spec_list, 0)[0]
        spec_list = [('distribcell', [0, 1, (6,1,2,1), 5, 2, (4,1,2,1), 1, 101])]
        outstr += '%s\n' % sp.get_value(1, spec_list, 0)[0]
        spec_list = [('distribcell', [0, 1, (6,2,1,1), 5, 2, (4,1,2,1), 1, 101])]
        outstr += '%s\n' % sp.get_value(1, spec_list, 0)[0]
        spec_list = [('distribcell', [0, 1, (6,2,2,1), 5, 2, (4,1,2,1), 1, 101])]
        outstr += '%s\n' % sp.get_value(1, spec_list, 0)[0]
        
        spec_list = [('distribcell', [0, 1, (6,1,1,1), 5, 2, (4,1,1,1), 2, 201])]
        outstr += '%s\n' % sp.get_value(2, spec_list, 0)[0]
        spec_list = [('distribcell', [0, 1, (6,1,1,1), 5, 2, (4,2,2,1), 2, 201])]
        outstr += '%s\n' % sp.get_value(2, spec_list, 0)[0]
        spec_list = [('distribcell', [0, 1, (6,2,1,1), 5, 2, (4,1,1,1), 2, 201])]
        outstr += '%s\n' % sp.get_value(2, spec_list, 0)[0]
        spec_list = [('distribcell', [0, 1, (6,2,1,1), 5, 2, (4,2,2,1), 2, 201])]
        outstr += '%s\n' % sp.get_value(2, spec_list, 0)[0]
        spec_list = [('distribcell', [0, 1, (6,1,2,1), 5, 2, (4,1,1,1), 2, 201])]
        outstr += '%s\n' % sp.get_value(2, spec_list, 0)[0]
        spec_list = [('distribcell', [0, 1, (6,1,2,1), 5, 2, (4,2,2,1), 2, 201])]
        outstr += '%s\n' % sp.get_value(2, spec_list, 0)[0]
        spec_list = [('distribcell', [0, 1, (6,2,2,1), 5, 2, (4,1,1,1), 2, 201])]
        outstr += '%s\n' % sp.get_value(2, spec_list, 0)[0]
        spec_list = [('distribcell', [0, 1, (6,2,2,1), 5, 2, (4,2,2,1), 2, 201])]
        outstr += '%s\n' % sp.get_value(2, spec_list, 0)[0]
        
        spec_list = [('distribcell', [0, 1, (6,1,1,1), 5, 2, (4,2,1,1), 3, 301])]
        outstr += '%s\n' % sp.get_value(3, spec_list, 0)[0]
        spec_list = [('distribcell', [0, 1, (6,1,2,1), 5, 2, (4,2,1,1), 3, 301])]
        outstr += '%s\n' % sp.get_value(3, spec_list, 0)[0]
        spec_list = [('distribcell', [0, 1, (6,2,1,1), 5, 2, (4,2,1,1), 3, 301])]
        outstr += '%s\n' % sp.get_value(3, spec_list, 0)[0]
        spec_list = [('distribcell', [0, 1, (6,2,2,1), 5, 2, (4,2,1,1), 3, 301])]
        outstr += '%s\n' % sp.get_value(3, spec_list, 0)[0]
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