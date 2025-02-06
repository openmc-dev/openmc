import os

import openmc
from openmc.examples import random_ray_three_region_cube

from tests.testing_harness import WeightWindowPyAPITestHarness


class MGXSTestHarness(WeightWindowPyAPITestHarness):
    def _cleanup(self):
        super()._cleanup()
        f = 'mgxs.h5'
        if os.path.exists(f):
            os.remove(f)


def test_random_ray_adjoint_fixed_source():
    model = random_ray_three_region_cube()

    # The base model has a resolution of 12, so we overlay
    # something else for FW-CADIS
    n = 15
    width = 30.0
    ww_mesh = openmc.RegularMesh()
    ww_mesh.dimension = (n, n, n)
    ww_mesh.lower_left = (0.0, 0.0, 0.0)
    ww_mesh.upper_right = (width, width, width)

    wwg = openmc.WeightWindowGenerator(
        method="fw_cadis", mesh=ww_mesh, max_realizations=model.settings.batches)
    model.settings.weight_window_generators = wwg
    model.settings.random_ray['volume_estimator'] = 'naive'
    
    root = model.geometry.root_universe
    model.settings.random_ray['source_region_meshes'] = [(ww_mesh, [root])]

    model.settings.particles = 400

    harness = MGXSTestHarness('statepoint.10.h5', model)
    harness.main()
