import os

import openmc
from openmc.examples import random_ray_three_region_cube

from tests.testing_harness import TolerantPyAPITestHarness


class MGXSTestHarness(TolerantPyAPITestHarness):
    def _cleanup(self):
        super()._cleanup()
        f = 'mgxs.h5'
        if os.path.exists(f):
            os.remove(f)


def test_random_ray_linear_source_stability():
    # A linear source run in the regime that stresses the batch-consistent
    # flux update: the naive volume estimator, an overlay source-region mesh
    # whose cells receive modest per-batch hit counts, and the example's
    # optically thin scattering-dominated interior.
    openmc.reset_auto_ids()
    model = random_ray_three_region_cube()
    model.settings.random_ray['source_shape'] = 'linear'
    model.settings.random_ray['volume_estimator'] = 'naive'
    mesh = openmc.RegularMesh()
    mesh.lower_left = (0.0, 0.0, 0.0)
    mesh.upper_right = (30.0, 30.0, 30.0)
    mesh.dimension = (12, 12, 12)
    model.settings.random_ray['source_region_meshes'] = [
        (mesh, [model.geometry.root_universe])]
    model.settings.inactive = 30
    model.settings.batches = 60
    harness = MGXSTestHarness('statepoint.60.h5', model)
    harness.main()
