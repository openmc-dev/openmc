import os
import glob

import openmc

from tests.testing_harness import PyAPITestHarness


class StatepointLatestPyAPITestHarness(PyAPITestHarness):
    def _test_output_created(self):
        """Make sure running statepoint files have been created/rotated."""
        # Call parent checks (writes inputs, etc.)
        super()._test_output_created()
        # Expect two running files (batches 3 and 4 when -2 is specified)
        running = glob.glob(os.path.join(os.getcwd(), 'statepoint.running.*.h5'))
        assert len(running) == 2, 'Expected 2 running statepoint files to exist.'


def test_statepoint_latest_pyapi():
    # Build a minimal model programmatically
    model = openmc.Model()

    mat = openmc.Material()
    mat.set_density('g/cm3', 1.0)
    mat.add_nuclide('H1', 1.0)
    model.materials.append(mat)

    # Simple box cell
    box = openmc.model.RectangularParallelepiped(-10, 10, -10, 10, -10, 10)
    cell = openmc.Cell(fill=mat, region=box)
    model.geometry = openmc.Geometry([cell])

    # Settings: small run to exercise feature
    model.settings.batches = 4
    model.settings.inactive = 0
    model.settings.particles = 100

    # Use Python API to keep running statepoints for last 2 batches (using negative batch number)
    model.settings.statepoint = {'batches': -2}

    harness = StatepointLatestPyAPITestHarness('statepoint.4.h5', model)
    harness.main()
