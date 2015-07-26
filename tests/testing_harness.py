from __future__ import print_function

import filecmp
import glob
import hashlib
from optparse import OptionParser
import os
import shutil
from subprocess import Popen, STDOUT, PIPE, call
import sys

import numpy as np

sys.path.insert(0, '../..')
from openmc.statepoint import StatePoint
from openmc.executor import Executor
import openmc.particle_restart as pr


class TestHarness(object):
    """General class for running OpenMC regression tests."""
    def __init__(self, statepoint_name, tallies_present=False):
        self._sp_name = statepoint_name
        self._tallies = tallies_present
        self._opts = None
        self._args = None

    def main(self):
        """Accept commandline arguments and either run or update tests."""
        self._parse_args()
        if self._opts.update:
            self.update_results()
        else:
            self.execute_test()

    def execute_test(self):
        """Run OpenMC with the appropriate arguments and check the outputs."""
        try:
            self._run_openmc()
            self._test_output_created()
            results = self._get_results()
            self._write_results(results)
            self._compare_results()
        finally:
            self._cleanup()

    def update_results(self):
        """Update the results_true using the current version of OpenMC."""
        try:
            self._run_openmc()
            self._test_output_created()
            results = self._get_results()
            self._write_results(results)
            self._overwrite_results()
        finally:
            self._cleanup()

    def _parse_args(self):
        parser = OptionParser()
        parser.add_option('--exe', dest='exe', default='openmc')
        parser.add_option('--mpi_exec', dest='mpi_exec', default=None)
        parser.add_option('--mpi_np', dest='mpi_np', type=int, default=3)
        parser.add_option('--update', dest='update', action='store_true',
                          default=False)
        (self._opts, self._args) = parser.parse_args()

    def _run_openmc(self):
        executor = Executor()

        if self._opts.mpi_exec is not None:
            returncode = executor.run_simulation(mpi_procs=self._opts.mpi_np,
                                                 openmc_exec=self._opts.exe,
                                                 mpi_exec=self._opts.mpi_exec)

        else:
            returncode = executor.run_simulation(openmc_exec=self._opts.exe)

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

    def _get_results(self, hash_output=False):
        """Digest info in the statepoint and return as a string."""
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

        # Hash the results if necessary.
        if hash_output:
            sha512 = hashlib.sha512()
            sha512.update(outstr.encode('utf-8'))
            outstr = sha512.hexdigest()

        return outstr

    def _write_results(self, results_string):
        """Write the results to an ASCII file."""
        with open('results_test.dat','w') as fh:
            fh.write(results_string)

    def _overwrite_results(self):
        """Overwrite the results_true with the results_test."""
        shutil.copyfile('results_test.dat', 'results_true.dat')

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


class HashedTestHarness(TestHarness):
    """Specialized TestHarness that hashes the results."""
    def _get_results(self):
        """Digest info in the statepoint and return as a string."""
        return TestHarness._get_results(self, True)


class PlotTestHarness(TestHarness):
    """Specialized TestHarness for running OpenMC plotting tests.""" 
    def __init__(self, plot_names):
        self._plot_names = plot_names
        self._opts = None
        self._args = None

    def _run_openmc(self):
        executor = Executor()
        returncode = executor.plot_geometry(openmc_exec=self._opts.exe)
        assert returncode == 0, 'OpenMC did not exit successfully.'

    def _test_output_created(self):
        """Make sure *.ppm has been created."""
        for fname in self._plot_names:
            assert os.path.exists(os.path.join(os.getcwd(), fname)), \
                 'Plot output file does not exist.'

    def _cleanup(self):
        TestHarness._cleanup(self)
        output = glob.glob(os.path.join(os.getcwd(), '*.ppm'))
        for f in output:
            if os.path.exists(f):
                os.remove(f)

    def _get_results(self):
        """Return a string hash of the plot files."""
        # Find the plot files.
        plot_files = glob.glob(os.path.join(os.getcwd(), '*.ppm'))

        # Read the plot files.
        outstr = bytes()
        for fname in sorted(plot_files):
            with open(fname, 'rb') as fh:
                outstr += fh.read()

        # Hash the information and return.
        sha512 = hashlib.sha512()
        sha512.update(outstr)
        outstr = sha512.hexdigest()

        return outstr


class CMFDTestHarness(TestHarness):
    """Specialized TestHarness for running OpenMC CMFD tests.""" 
    def _get_results(self):
        """Digest info in the statepoint and return as a string."""
        # Read the statepoint file.
        statepoint = glob.glob(os.path.join(os.getcwd(), self._sp_name))[0]
        sp = StatePoint(statepoint)
        sp.read_results()

        # Write out the eigenvalue and tallies.
        outstr = TestHarness._get_results(self)

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

        return outstr


class ParticleRestartTestHarness(TestHarness):
    """Specialized TestHarness for running OpenMC particle restart tests.""" 
    def _test_output_created(self):
        """Make sure the restart file has been created."""
        particle = glob.glob(os.path.join(os.getcwd(), self._sp_name))
        assert len(particle) == 1, 'Either multiple or no particle restart ' \
             'files exist.'
        assert particle[0].endswith('binary') \
             or particle[0].endswith('h5'), \
             'Particle restart file is not a binary or hdf5 file.'

    def _get_results(self):
        """Digest info in the statepoint and return as a string."""
        # Read the particle restart file.
        particle = glob.glob(os.path.join(os.getcwd(), self._sp_name))[0]
        p = pr.Particle(particle)

        # Write out the properties.
        outstr = ''
        outstr += 'current batch:\n'
        outstr += "{0:12.6E}\n".format(p.current_batch)
        outstr += 'current gen:\n'
        outstr += "{0:12.6E}\n".format(p.current_gen)
        outstr += 'particle id:\n'
        outstr += "{0:12.6E}\n".format(p.id)
        outstr += 'run mode:\n'
        outstr += "{0:12.6E}\n".format(p.run_mode)
        outstr += 'particle weight:\n'
        outstr += "{0:12.6E}\n".format(p.weight)
        outstr += 'particle energy:\n'
        outstr += "{0:12.6E}\n".format(p.energy)
        outstr += 'particle xyz:\n'
        outstr += "{0:12.6E} {1:12.6E} {2:12.6E}\n".format(p.xyz[0],p.xyz[1],
                                                           p.xyz[2])
        outstr += 'particle uvw:\n'
        outstr += "{0:12.6E} {1:12.6E} {2:12.6E}\n".format(p.uvw[0],p.uvw[1],
                                                           p.uvw[2])

        return outstr
