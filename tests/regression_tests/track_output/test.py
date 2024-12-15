import glob
import os
from pathlib import Path
from subprocess import call

import numpy as np
import openmc
import pytest

from tests.testing_harness import TestHarness, config


class TrackTestHarness(TestHarness):
    def _test_output_created(self):
        """Make sure statepoint.* and track* have been created."""
        TestHarness._test_output_created(self)

        if config["mpi"] and int(config["mpi_np"]) > 1:
            outputs = Path.cwd().glob("tracks_p*.h5")
            assert len(list(outputs)) == int(config["mpi_np"])
        else:
            assert Path("tracks.h5").is_file()

    def _get_results(self):
        """Get data from track file and return as a string."""

        # For MPI mode, combine track files
        if config["mpi"]:
            call(
                ["../../../scripts/openmc-track-combine", "-o", "tracks.h5"]
                + glob.glob("tracks_p*.h5")
            )

        # Get string of track file information
        outstr = ""
        tracks = openmc.Tracks("tracks.h5")
        for track in tracks:
            with np.printoptions(formatter={"float_kind": "{:.6e}".format}):
                for ptrack in track:
                    outstr += f"{ptrack.particle} {ptrack.states}\n"

        return outstr

    def _cleanup(self):
        TestHarness._cleanup(self)
        output = glob.glob("tracks*") + glob.glob("poly*")
        for f in output:
            if os.path.exists(f):
                os.remove(f)


def test_track_output():
    # If vtk python module is not available, we can't run track.py so skip this
    # test.
    vtk = pytest.importorskip("vtk")
    harness = TrackTestHarness("statepoint.2.h5")
    harness.main()
