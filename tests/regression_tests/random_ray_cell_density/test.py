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


@pytest.mark.parametrize("run_mode", ["eigen", "fs"])
def test_random_ray_basic(run_mode):
    with change_directory(run_mode):
        if run_mode == "eigen":
            openmc.reset_auto_ids()
            model = random_ray_lattice()
            # Double the densities of the lower-left fuel pin -> cell instances [0, 8).
            for id, cell in model.geometry.get_all_cells().items():
                if cell.fill.name == "UO2 fuel":
                    cell.density = [((i < 8) + 1.0) for i in range(24)]

            # Gold file was generated with manually scaled fuel cross sections.
            harness = MGXSTestHarness('statepoint.10.h5', model)
            harness.main()
        else:
            openmc.reset_auto_ids()
            model = random_ray_three_region_cube()
            # Increase the density in the source region.
            for id, cell in model.geometry.get_all_cells().items():
                if cell.fill.name == "source":
                    cell.density = 1e3

            # Gold file was generated with manually scaled source cross sections.
            harness = MGXSTestHarness('statepoint.10.h5', model)
            harness.main()
