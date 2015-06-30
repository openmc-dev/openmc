#!/usr/bin/env python

import shutil
import sys

sys.path.insert(0, '..')
from testing_harness import *


class TrackTestHarness(TestHarness):
    def _test_output_created(self):
        """Make sure statepoint.* and track* have been created."""
        TestHarness._test_output_created(self)

        outputs = [glob.glob(''.join((os.getcwd(), '/track_1_1_1.*')))]
        outputs.append(glob.glob(''.join((os.getcwd(), '/track_1_1_2.*'))))
        for files in outputs:
            assert len(files) == 1, 'Multiple or no track files detected.'
            assert files[0].endswith('binary') or files[0].endswith('h5'),\
            'Track files are not binary or hdf5 files'

    def _get_results(self):
        """Digest info in the statepoint and return as a string."""
        # Run the track-to-vtk conversion script.
        call(['../../scripts/openmc-track-to-vtk', '-o', 'poly'] +
             glob.glob(''.join((os.getcwd(), '/track*'))))

        # Make sure the vtk file was created then return it's contents.
        poly = os.path.join(os.getcwd(), 'poly.pvtp')
        assert os.path.isfile(poly), 'poly.pvtp file not found.'

        with open(poly) as fin:
            outstr = fin.read()

        return outstr

    def _cleanup(self):
        TestHarness._cleanup(self)
        output = glob.glob(os.path.join(os.getcwd(), 'track*'))
        output += glob.glob(os.path.join(os.getcwd(), 'poly*'))
        for f in output:
            if os.path.exists(f):
                os.remove(f)


if __name__ == '__main__':
    # If vtk python module is not available, we can't run track.py so skip this
    # test.
    try:
        import vtk
    except ImportError:
        print('----------------Skipping test-------------')
        shutil.copy('results_true.dat', 'results_test.dat')
        exit()
    harness = TrackTestHarness('statepoint.2.*')
    harness.main()
