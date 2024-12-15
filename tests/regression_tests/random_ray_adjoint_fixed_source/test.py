import os

from openmc.examples import random_ray_three_region_cube

from tests.testing_harness import TolerantPyAPITestHarness


class MGXSTestHarness(TolerantPyAPITestHarness):
    def _cleanup(self):
        super()._cleanup()
        f = "mgxs.h5"
        if os.path.exists(f):
            os.remove(f)


def test_random_ray_adjoint_fixed_source():
    model = random_ray_three_region_cube()
    model.settings.random_ray["adjoint"] = True
    harness = MGXSTestHarness("statepoint.10.h5", model)
    harness.main()
