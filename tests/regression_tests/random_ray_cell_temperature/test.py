import os

import openmc
from openmc.examples import random_ray_lattice, random_ray_three_region_cube
from openmc.utility_funcs import change_directory
import pytest

from tests.testing_harness import TolerantPyAPITestHarness


class MGXSTestHarness(TolerantPyAPITestHarness):
    def _cleanup(self):
        super()._cleanup()
        f = 'mgxs.h5'
        if os.path.exists(f):
            os.remove(f)


def test_random_ray_basic():
    openmc.reset_auto_ids()
    model = random_ray_lattice(second_temp=True)
    # Set the temperature of the lower-left pin to 395 K -> cell instances [0, 8).
    # All other pins are set to 295.
    for id, cell in model.geometry.get_all_cells().items():
        if cell.fill.name == "UO2 fuel":
            cell.temperature = [(100.0 * (i < 8) + 295.0) for i in range(24)]

    model.settings.temperature = {
        'method' : 'nearest',
        'tolerance' : 10.0,
        'range' : (200.0, 400.0)
    }

    # Gold file was generated with manually scaled fuel cross sections.
    harness = MGXSTestHarness('statepoint.10.h5', model)
    harness.main()
