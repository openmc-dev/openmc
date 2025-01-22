from difflib import unified_diff
from subprocess import check_call
import filecmp
import glob
import h5py
import hashlib
import os
import shutil

import numpy as np
import openmc
from openmc.examples import pwr_core
from colorama import Fore, init

from tests.regression_tests import config

init()


def colorize(diff):
    """Produce colored diff for test results"""
    for line in diff:
        if line.startswith('+'):
            yield Fore.RED + line + Fore.RESET
        elif line.startswith('-'):
            yield Fore.GREEN + line + Fore.RESET
        elif line.startswith('^'):
            yield Fore.BLUE + line + Fore.RESET
        else:
            yield line


class TestHarness:
    """General class for running OpenMC regression tests."""

    def __init__(self, statepoint_name):
        self._sp_name = statepoint_name

    def main(self):
        """Accept commandline arguments and either run or update tests."""
        if config['update']:
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

    def _run_openmc(self):
        if config['mpi']:
            mpi_args = [config['mpiexec'], '-n', config['mpi_np']]
            openmc.run(openmc_exec=config['exe'], mpi_args=mpi_args,
              event_based=config['event'])
        else:
            openmc.run(openmc_exec=config['exe'], event_based=config['event'])

    def _test_output_created(self):
        """Make sure statepoint.* and tallies.out have been created."""
        statepoint = glob.glob(self._sp_name)
        assert len(statepoint) == 1, 'Either multiple or no statepoint files' \
            ' exist.'
        assert statepoint[0].endswith('h5'), \
            'Statepoint file is not a HDF5 file.'
        if os.path.exists('tallies.xml'):
            assert os.path.exists('tallies.out'), \
                'Tally output file does not exist.'

    def _get_results(self, hash_output=False):
        """Digest info in the statepoint and return as a string."""
        # Read the statepoint file.
        statepoint = glob.glob(self._sp_name)[0]
        with openmc.StatePoint(statepoint) as sp:
            outstr = ''
            if sp.run_mode == 'eigenvalue':
                # Write out k-combined.
                outstr += 'k-combined:\n'
                form = '{0:12.6E} {1:12.6E}\n'
                outstr += form.format(sp.keff.n, sp.keff.s)

            # Write out tally data.
            for i, tally_ind in enumerate(sp.tallies):
                tally = sp.tallies[tally_ind]
                results = np.zeros((tally.sum.size * 2, ))
                results[0::2] = tally.sum.ravel()
                results[1::2] = tally.sum_sq.ravel()
                results = ['{0:12.6E}'.format(x) for x in results]

                outstr += 'tally {}:\n'.format(i + 1)
                outstr += '\n'.join(results) + '\n'

        # Hash the results if necessary.
        if hash_output:
            sha512 = hashlib.sha512()
            sha512.update(outstr.encode('utf-8'))
            outstr = sha512.hexdigest()

        return outstr

    @property
    def statepoint_name(self):
        return self._sp_name

    def _write_results(self, results_string):
        """Write the results to an ASCII file."""
        with open('results_test.dat', 'w') as fh:
            fh.write(results_string)

    def _overwrite_results(self):
        """Overwrite the results_true with the results_test."""
        shutil.copyfile('results_test.dat', 'results_true.dat')

    def _compare_results(self):
        """Make sure the current results agree with the reference."""
        compare = filecmp.cmp('results_test.dat', 'results_true.dat')
        if not compare:
            expected = open('results_true.dat').readlines()
            actual = open('results_test.dat').readlines()
            diff = unified_diff(expected, actual, 'results_true.dat',
                                'results_test.dat')
            print('Result differences:')
            print(''.join(colorize(diff)))
            os.rename('results_test.dat', 'results_error.dat')
        assert compare, 'Results do not agree'

    def _cleanup(self):
        """Delete statepoints, tally, and test files."""
        output = glob.glob('statepoint.*.h5')
        output += ['tallies.out', 'results_test.dat', 'summary.h5']
        output += glob.glob('volume_*.h5')
        for f in output:
            if os.path.exists(f):
                os.remove(f)


class HashedTestHarness(TestHarness):
    """Specialized TestHarness that hashes the results."""

    def _get_results(self):
        """Digest info in the statepoint and return as a string."""
        return super()._get_results(True)


class CMFDTestHarness(TestHarness):
    """Specialized TestHarness for running OpenMC CMFD tests."""

    def __init__(self, statepoint_name, cmfd_run):
        super().__init__(statepoint_name)
        self._create_cmfd_result_str(cmfd_run)

    def _create_cmfd_result_str(self, cmfd_run):
        """Create CMFD result string from variables of CMFDRun instance"""
        outstr = 'cmfd indices\n'
        outstr += '\n'.join(['{:.6E}'.format(x) for x in cmfd_run.indices])
        outstr += '\nk cmfd\n'
        outstr += '\n'.join(['{:.6E}'.format(x) for x in cmfd_run.k_cmfd])
        outstr += '\ncmfd entropy\n'
        outstr += '\n'.join(['{:.6E}'.format(x) for x in cmfd_run.entropy])
        outstr += '\ncmfd balance\n'
        outstr += '\n'.join(['{:.5E}'.format(x) for x in cmfd_run.balance])
        outstr += '\ncmfd dominance ratio\n'
        outstr += '\n'.join(['{:.3E}'.format(x) for x in cmfd_run.dom])
        outstr += '\ncmfd openmc source comparison\n'
        outstr += '\n'.join(['{:.6E}'.format(x) for x in cmfd_run.src_cmp])
        outstr += '\ncmfd source\n'
        cmfdsrc = np.reshape(cmfd_run.cmfd_src, np.prod(cmfd_run.indices),
                             order='F')
        outstr += '\n'.join(['{:.6E}'.format(x) for x in cmfdsrc])
        outstr += '\n'
        self._cmfdrun_results = outstr

    def execute_test(self):
        """Don't call _run_openmc as OpenMC will be called through C API for
        CMFD tests, and write CMFD results that were passsed as argument

        """
        try:
            self._test_output_created()
            results = self._get_results()
            results += self._cmfdrun_results
            self._write_results(results)
            self._compare_results()
        finally:
            self._cleanup()

    def update_results(self):
        """Don't call _run_openmc as OpenMC will be called through C API for
        CMFD tests, and write CMFD results that were passsed as argument

        """
        try:
            self._test_output_created()
            results = self._get_results()
            results += self._cmfdrun_results
            self._write_results(results)
            self._overwrite_results()
        finally:
            self._cleanup()

    def _cleanup(self):
        """Delete output files for numpy matrices and flux vectors."""
        super()._cleanup()
        output = ['loss.npz', 'loss.dat', 'prod.npz', 'prod.dat',
                  'fluxvec.npy', 'fluxvec.dat']
        for f in output:
            if os.path.exists(f):
                os.remove(f)


class ParticleRestartTestHarness(TestHarness):
    """Specialized TestHarness for running OpenMC particle restart tests."""

    def _run_openmc(self):
        # Set arguments
        args = {'openmc_exec': config['exe']}
        if config['mpi']:
            args['mpi_args'] = [config['mpiexec'], '-n', config['mpi_np']]

        # Initial run
        openmc.run(**args)

        # Run particle restart
        args.update({'restart_file': self._sp_name})
        openmc.run(**args)

    def _test_output_created(self):
        """Make sure the restart file has been created."""
        particle = glob.glob(self._sp_name)
        assert len(particle) == 1, 'Either multiple or no particle restart ' \
            'files exist.'
        assert particle[0].endswith('h5'), \
            'Particle restart file is not a HDF5 file.'

    def _get_results(self):
        """Digest info in the statepoint and return as a string."""
        # Read the particle restart file.
        particle = glob.glob(self._sp_name)[0]
        p = openmc.Particle(particle)

        # Write out the properties.
        outstr = ''
        outstr += 'current batch:\n'
        outstr += "{0:12.6E}\n".format(p.current_batch)
        outstr += 'current generation:\n'
        outstr += "{0:12.6E}\n".format(p.current_generation)
        outstr += 'particle id:\n'
        outstr += "{0:12.6E}\n".format(p.id)
        outstr += 'run mode:\n'
        outstr += "{0}\n".format(p.run_mode)
        outstr += 'particle weight:\n'
        outstr += "{0:12.6E}\n".format(p.weight)
        outstr += 'particle energy:\n'
        outstr += "{0:12.6E}\n".format(p.energy)
        outstr += 'particle xyz:\n'
        outstr += "{0:12.6E} {1:12.6E} {2:12.6E}\n".format(p.xyz[0], p.xyz[1],
                                                           p.xyz[2])
        outstr += 'particle uvw:\n'
        outstr += "{0:12.6E} {1:12.6E} {2:12.6E}\n".format(p.uvw[0], p.uvw[1],
                                                           p.uvw[2])

        return outstr

    def _cleanup(self):
        """Delete particle restart files."""
        super()._cleanup()
        output = glob.glob('particle*.h5')
        for f in output:
            os.remove(f)


class PyAPITestHarness(TestHarness):
    def __init__(self, statepoint_name, model=None, inputs_true=None):
        super().__init__(statepoint_name)
        if model is None:
            self._model = pwr_core()
        else:
            self._model = model
        self._model.plots = []

        self.inputs_true = "inputs_true.dat" if not inputs_true else inputs_true

    def main(self):
        """Accept commandline arguments and either run or update tests."""
        if config['build_inputs']:
            self._build_inputs()
        elif config['update']:
            self.update_results()
        else:
            self.execute_test()

    def execute_test(self):
        """Build input XMLs, run OpenMC, and verify correct results."""
        try:
            self._build_inputs()
            inputs = self._get_inputs()
            self._write_inputs(inputs)
            self._compare_inputs()
            self._run_openmc()
            self._test_output_created()
            results = self._get_results()
            self._write_results(results)
            self._compare_results()
        finally:
            self._cleanup()

    def update_results(self):
        """Update results_true.dat and inputs_true.dat"""
        try:
            self._build_inputs()
            inputs = self._get_inputs()
            self._write_inputs(inputs)
            self._overwrite_inputs()
            self._run_openmc()
            self._test_output_created()
            results = self._get_results()
            self._write_results(results)
            self._overwrite_results()
        finally:
            self._cleanup()

    def _build_inputs(self):
        """Write input XML files."""
        self._model.export_to_model_xml()

    def _get_inputs(self):
        """Return a hash digest of the input XML files."""
        xmls = ['model.xml', 'plots.xml']
        return ''.join([open(fname).read() for fname in xmls
                        if os.path.exists(fname)])

    def _write_inputs(self, input_digest):
        """Write the digest of the input XMLs to an ASCII file."""
        with open('inputs_test.dat', 'w') as fh:
            fh.write(input_digest)

    def _overwrite_inputs(self):
        """Overwrite inputs_true.dat with inputs_test.dat"""
        shutil.copyfile('inputs_test.dat', self.inputs_true)

    def _compare_inputs(self):
        """Make sure the current inputs agree with the _true standard."""
        compare = filecmp.cmp('inputs_test.dat', self.inputs_true)
        if not compare:
            expected = open(self.inputs_true, 'r').readlines()
            actual = open('inputs_test.dat', 'r').readlines()
            diff = unified_diff(expected, actual, self.inputs_true,
                                'inputs_test.dat')
            print('Input differences:')
            print(''.join(colorize(diff)))
            os.rename('inputs_test.dat', 'inputs_error.dat')
        assert compare, 'Input files are broken.'

    def _cleanup(self):
        """Delete XMLs, statepoints, tally, and test files."""
        super()._cleanup()
        output = ['materials.xml', 'geometry.xml', 'settings.xml',
                  'tallies.xml', 'plots.xml', 'inputs_test.dat', 'model.xml']
        for f in output:
            if os.path.exists(f):
                os.remove(f)


class HashedPyAPITestHarness(PyAPITestHarness):
    def _get_results(self):
        """Digest info in the statepoint and return as a string."""
        return super()._get_results(True)


class TolerantPyAPITestHarness(PyAPITestHarness):
    """Specialized harness for running tests that involve significant levels
    of floating point non-associativity when using shared memory parallelism
    due to single precision usage (e.g., as in the random ray solver).

    """
    def _are_files_equal(self, actual_path, expected_path, tolerance):
        def isfloat(value):
            try:
                float(value)
                return True
            except ValueError:
                return False

        def tokenize(line):
            return line.strip().split()

        def compare_tokens(token1, token2):
            if isfloat(token1) and isfloat(token2):
                float1, float2 = float(token1), float(token2)
                return abs(float1 - float2) <= tolerance * max(abs(float1), abs(float2))
            else:
                return token1 == token2

        expected = open(expected_path).readlines()
        actual = open(actual_path).readlines()

        if len(expected) != len(actual):
            return False

        for line1, line2 in zip(expected, actual):
            tokens1 = tokenize(line1)
            tokens2 = tokenize(line2)

            if len(tokens1) != len(tokens2):
                return False

            for token1, token2 in zip(tokens1, tokens2):
                if not compare_tokens(token1, token2):
                    return False

        return True

    def _compare_results(self):
        """Make sure the current results agree with the reference."""
        compare = self._are_files_equal('results_test.dat', 'results_true.dat', 1e-6)
        if not compare:
            expected = open('results_true.dat').readlines()
            actual = open('results_test.dat').readlines()
            diff = unified_diff(expected, actual, 'results_true.dat',
                                'results_test.dat')
            print('Result differences:')
            print(''.join(colorize(diff)))
            os.rename('results_test.dat', 'results_error.dat')
        assert compare, 'Results do not agree'


class WeightWindowPyAPITestHarness(PyAPITestHarness):
    def _get_results(self):
        """Digest info in the weight window file and return as a string."""
        ww = openmc.hdf5_to_wws()[0]

        # Access the weight window bounds
        lower_bound = ww.lower_ww_bounds
        upper_bound = ww.upper_ww_bounds

        # Flatten both arrays
        flattened_lower_bound = lower_bound.flatten()
        flattened_upper_bound = upper_bound.flatten()

        # Convert each element to a string in scientific notation with 2 decimal places
        formatted_lower_bound = [f'{x:.2e}' for x in flattened_lower_bound]
        formatted_upper_bound = [f'{x:.2e}' for x in flattened_upper_bound]

        # Concatenate the formatted arrays
        concatenated_strings = [
            "Lower Bounds"] + formatted_lower_bound + ["Upper Bounds"] + formatted_upper_bound

        # Join the concatenated strings into a single string with newline characters
        final_string = '\n'.join(concatenated_strings)

        # Prepend the mesh text description and return final string
        final_string = str(ww.mesh) + final_string
        return final_string
    
    def _cleanup(self):
        super()._cleanup()
        f = 'weight_windows.h5'
        if os.path.exists(f):
            os.remove(f)


class PlotTestHarness(TestHarness):
    """Specialized TestHarness for running OpenMC plotting tests."""
    def __init__(self, plot_names, voxel_convert_checks=[]):
        super().__init__(None)
        self._plot_names = plot_names
        self._voxel_convert_checks = voxel_convert_checks

    def _run_openmc(self):
        openmc.plot_geometry(openmc_exec=config['exe'])

        # Check that voxel h5 can be converted to vtk
        for voxel_h5_filename in self._voxel_convert_checks:
            check_call(['../../../scripts/openmc-voxel-to-vtk'] +
                       glob.glob(voxel_h5_filename))

    def _test_output_created(self):
        """Make sure *.png has been created."""
        for fname in self._plot_names:
            assert os.path.exists(fname), 'Plot output file does not exist.'

    def _cleanup(self):
        super()._cleanup()
        for fname in self._plot_names:
            if os.path.exists(fname):
                os.remove(fname)

    def _get_results(self):
        """Return a string hash of the plot files."""
        outstr = bytes()

        for fname in self._plot_names:
            if fname.endswith('.png'):
                # Add PNG output to results
                with open(fname, 'rb') as fh:
                    outstr += fh.read()
            elif fname.endswith('.h5'):
                # Add voxel data to results
                with h5py.File(fname, 'r') as fh:
                    outstr += fh.attrs['filetype']
                    outstr += fh.attrs['num_voxels'].tobytes()
                    outstr += fh.attrs['lower_left'].tobytes()
                    outstr += fh.attrs['voxel_width'].tobytes()
                    outstr += fh['data'][()].tobytes()

        # Hash the information and return.
        sha512 = hashlib.sha512()
        sha512.update(outstr)
        outstr = sha512.hexdigest()

        return outstr
