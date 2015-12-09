#!/usr/bin/env python

import os
import sys
import glob

sys.path.insert(0, os.pardir)
from testing_harness import TestHarness
from openmc.executor import Executor
import h5py
import shutil

class DomainDecompOTFMatsTestHarness(TestHarness):
    def main(self):
        """Accept commandline arguments and either run or update tests."""
        self._parse_args()
        
        if self._opts.mpi_exec is None:
            # This test is only relevant if MPI is enabled
            return
        shutil.copy('materials.h5.keep','materials.h5')
        if self._opts.update:
            self.update_results()
        else:
            self.execute_test()
    
    def _run_openmc(self):
        executor = Executor()
        returncode = executor.run_simulation(mpi_procs=4,
                                             openmc_exec=self._opts.exe,
                                             mpi_exec=self._opts.mpi_exec)
        assert returncode == 0, 'OpenMC did not exit successfully.'
    
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
    harness = DomainDecompOTFMatsTestHarness('statepoint.20.*')
    harness.main()

