import glob
import hashlib
import os
from subprocess import check_call

import h5py
import openmc
import pytest

from tests.testing_harness import TestHarness
from tests.regression_tests import config


vtk = pytest.importorskip('vtk')
class PlotVoxelTestHarness(TestHarness):
    """Specialized TestHarness for running OpenMC voxel plot tests."""
    def __init__(self, plot_names):
        super().__init__(None)
        self._plot_names = plot_names

    def _run_openmc(self):
        openmc.plot_geometry(openmc_exec=config['exe'])

        check_call(['../../../scripts/openmc-voxel-to-vtk'] +
                   glob.glob('plot_4.h5'))
                
    def _test_output_created(self):
        """Make sure *.ppm has been created."""
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
            if fname.endswith('.h5'):
                # Add voxel data to results
                with h5py.File(fname, 'r') as fh:
                    outstr += fh.attrs['filetype']
                    outstr += fh.attrs['num_voxels'].tostring()
                    outstr += fh.attrs['lower_left'].tostring()
                    outstr += fh.attrs['voxel_width'].tostring()
                    outstr += fh['data'].value.tostring()

        # Hash the information and return.
        sha512 = hashlib.sha512()
        sha512.update(outstr)
        outstr = sha512.hexdigest()

        return outstr


def test_plot_voxel():
    harness = PlotVoxelTestHarness(('plot_4.h5', 'plot.vti'))
    harness.main()
