import glob
import os
from subprocess import call
import shutil

import pytest

from tests.testing_harness import TestHarness


class TrackTestHarness(TestHarness):
    def _test_output_created(self):
        """Make sure statepoint.* and track* have been created."""
        TestHarness._test_output_created(self)

        outputs = glob.glob('track_1_1_*.h5')
        assert len(outputs) == 2, 'Expected two track files.'

    def _get_results(self):
        """Digest info in the statepoint and return as a string."""
        # Run the track-to-vtk conversion script.
        call(['../../../scripts/openmc-track-to-vtk', '-o', 'poly'] +
             glob.glob('track_1_1_*.h5'))

        # Make sure the vtk file was created then return it's contents.
        assert os.path.isfile('poly.pvtp'), 'poly.pvtp file not found.'

        with open('poly.pvtp', 'r') as fin:
            outstr = fin.read()

        return outstr

    def _cleanup(self):
        TestHarness._cleanup(self)
        output = glob.glob('track*') + glob.glob('poly*')
        for f in output:
            if os.path.exists(f):
                os.remove(f)


def test_track_output():
    # If vtk python module is not available, we can't run track.py so skip this
    # test.
    vtk = pytest.importorskip('vtk')
    harness = TrackTestHarness('statepoint.2.h5')
    harness.main()
