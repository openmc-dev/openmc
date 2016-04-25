#!/usr/bin/env python

import glob
import hashlib
import os
import sys
sys.path.insert(0, os.pardir)
from testing_harness import TestHarness

import h5py

import openmc


class PlotTestHarness(TestHarness):
    """Specialized TestHarness for running OpenMC plotting tests."""
    def __init__(self, plot_names):
        super(PlotTestHarness, self).__init__(None, False)
        self._plot_names = plot_names

    def _run_openmc(self):
        returncode = openmc.plot_geometry(openmc_exec=self._opts.exe)
        assert returncode == 0, 'OpenMC did not exit successfully.'

    def _test_output_created(self):
        """Make sure *.ppm has been created."""
        for fname in self._plot_names:
            assert os.path.exists(os.path.join(os.getcwd(), fname)), \
                 'Plot output file does not exist.'

    def _cleanup(self):
        super(PlotTestHarness, self)._cleanup()
        for fname in self._plot_names:
            path = os.path.join(os.getcwd(), fname)
            if os.path.exists(path):
                #os.remove(path)
                pass

    def _get_results(self):
        """Return a string hash of the plot files."""
        outstr = bytes()

        # Add PPM output to results
        ppm_files = glob.glob(os.path.join(os.getcwd(), '*.ppm'))
        for fname in sorted(ppm_files):
            with open(fname, 'rb') as fh:
                outstr += fh.read()

        # Add voxel data to results
        voxel_files = glob.glob(os.path.join(os.getcwd(), '*.voxel'))
        for fname in sorted(voxel_files):
            with h5py.File(fname, 'r') as fh:
                outstr += fh['filetype'].value
                outstr += fh['num_voxels'].value.tostring()
                outstr += fh['lower_left'].value.tostring()
                outstr += fh['voxel_width'].value.tostring()
                outstr += fh['data'].value.tostring()

        # Hash the information and return.
        sha512 = hashlib.sha512()
        sha512.update(outstr)
        outstr = sha512.hexdigest()

        return outstr


if __name__ == '__main__':
    harness = PlotTestHarness(('1_plot.ppm', '2_plot.ppm', '3_plot.ppm',
                               '4_plot.voxel'))
    harness.main()
