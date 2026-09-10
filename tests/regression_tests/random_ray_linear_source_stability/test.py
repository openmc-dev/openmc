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
    # flux update: the naive volume estimator updates every hit region with
    # its own batch volume, so the consistency term applies throughout, and
    # a ray count low enough to starve the example's 1728 regions (about 4%
    # of them missed per batch) makes the batch centroids scatter widely
    # about the accumulated ones.
    openmc.reset_auto_ids()
    model = random_ray_three_region_cube()
    model.settings.random_ray['source_shape'] = 'linear'
    model.settings.random_ray['volume_estimator'] = 'naive'
    model.settings.particles = 20
    model.settings.inactive = 30
    model.settings.batches = 60
    harness = MGXSTestHarness('statepoint.60.h5', model)
    harness.main()
