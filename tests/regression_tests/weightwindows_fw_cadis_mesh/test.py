import os

import openmc
from openmc.utility_funcs import change_directory
from openmc.examples import random_ray_three_region_cube
import pytest

from tests.testing_harness import WeightWindowPyAPITestHarness


class MGXSTestHarness(WeightWindowPyAPITestHarness):
    def _cleanup(self):
        super()._cleanup()
        f = 'mgxs.h5'
        if os.path.exists(f):
            os.remove(f)


@pytest.mark.parametrize("shape", ["flat", "linear"])
def test_weight_windows_fw_cadis_mesh(shape):
    with change_directory(shape):
        openmc.reset_auto_ids()

        model = random_ray_three_region_cube()

        # The base model has a resolution of 12, so we overlay
        # something else for FW-CADIS
        n = 14
        width = 30.0
        ww_mesh = openmc.RegularMesh()
        ww_mesh.dimension = (n, n, n)
        ww_mesh.lower_left = (0.0, 0.0, 0.0)
        ww_mesh.upper_right = (width, width, width)

        wwg = openmc.WeightWindowGenerator(
            method="fw_cadis", mesh=ww_mesh, max_realizations=model.settings.batches)
        model.settings.weight_window_generators = wwg
        
        root = model.geometry.root_universe
        model.settings.random_ray['source_region_meshes'] = [(ww_mesh, [root])]

        model.settings.particles = 150
        model.settings.batches = 14
        model.settings.inactive = 12

        model.settings.random_ray['source_shape'] = shape

        harness = MGXSTestHarness('statepoint.14.h5', model)
        harness.main()
