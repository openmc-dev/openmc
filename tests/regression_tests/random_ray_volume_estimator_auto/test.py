import os

import openmc
from openmc.utility_funcs import change_directory
from openmc.examples import random_ray_three_region_cube
import pytest

from tests.testing_harness import TolerantPyAPITestHarness


class MGXSTestHarness(TolerantPyAPITestHarness):
    def _cleanup(self):
        super()._cleanup()
        f = 'mgxs.h5'
        if os.path.exists(f):
            os.remove(f)


# The default "auto" volume estimator resolves by solve type: standard
# solves receive the adaptive estimator, while solves whose results feed
# variance reduction (any adjoint workflow, and weight window generation)
# receive the strict adaptive estimator. Neither case sets an estimator
# explicitly, so these golds pin the routing itself: if the resolution
# policy regresses, the affected case's results shift.
@pytest.mark.parametrize("solve", ["forward", "adjoint"])
def test_random_ray_volume_estimator_auto(solve):
    with change_directory(solve):
        openmc.reset_auto_ids()
        model = random_ray_three_region_cube()
        if solve == "adjoint":
            model.settings.random_ray['adjoint'] = True
        model.settings.particles = 10
        model.settings.inactive = 20
        model.settings.batches = 40

        harness = MGXSTestHarness('statepoint.40.h5', model)
        harness.main()
