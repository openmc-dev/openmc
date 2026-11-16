import os

import openmc
from openmc.examples import random_ray_lattice, random_ray_three_region_cube
from openmc.utility_funcs import change_directory
import pytest

from tests.testing_harness import VTKTestHarness


class MGXSTestHarness(VTKTestHarness):
    def _cleanup(self):
        super()._cleanup()
        f = 'mgxs.h5'
        if os.path.exists(f):
            os.remove(f)


@pytest.mark.parametrize("shape", ["flat", "linear"])
@pytest.mark.parametrize("run_mode", ["eigen", "fs"])
def test_random_ray_vtk(run_mode, shape):
    with change_directory(f'{run_mode}_{shape}'):
        openmc.reset_auto_ids()
        plot = openmc.VoxelPlot()
        if run_mode == 'eigen':
            model = random_ray_lattice()
            # Voxel centers off the azimuthal planes and rings of the pins
            plot.width = (2.52, 2.2, 1.0)
            plot.pixels = (4, 4, 1)
        else:
            model = random_ray_three_region_cube()
            # Voxel centers off the 2.5 cm lattice planes
            plot.origin = (12.75, 12.75, 12.75)
            plot.width = (24.9, 24.9, 24.9)
            plot.pixels = (3, 3, 3)
        model.settings.random_ray['source_shape'] = shape
        # The linear source keeps its gradients at zero for the first ten
        # batches, so put all the active batches after them
        model.settings.inactive = 10
        model.settings.batches = 15
        harness = MGXSTestHarness('statepoint.15.h5', model, plot)
        harness.main()
