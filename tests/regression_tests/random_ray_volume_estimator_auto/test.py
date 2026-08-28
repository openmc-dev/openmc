import os

import openmc
from openmc.utility_funcs import change_directory
from openmc.examples import random_ray_three_region_cube
import pytest

from tests.testing_harness import TolerantPyAPITestHarness


class MGXSTestHarness(TolerantPyAPITestHarness):
    def _cleanup(self):
        super()._cleanup()
        for f in ('mgxs.h5', 'weight_windows.h5'):
            if os.path.exists(f):
                os.remove(f)


# The default "auto" volume estimator resolves by solve type: standard
# solves receive the adaptive estimator, while solves whose results feed
# variance reduction (any adjoint workflow, and weight window generation)
# receive the strict adaptive estimator. No case sets an estimator
# explicitly, so these golds pin the routing itself: if the resolution
# policy regresses, the affected case's results shift. The forward and
# adjoint cases pin the adjoint-flag trigger; the weight_windows case pins
# the generator-presence trigger.
@pytest.mark.parametrize("solve", ["forward", "adjoint", "weight_windows"])
def test_random_ray_volume_estimator_auto(solve):
    with change_directory(solve):
        openmc.reset_auto_ids()
        model = random_ray_three_region_cube()
        if solve == "adjoint":
            model.settings.random_ray['adjoint'] = True
        elif solve == "weight_windows":
            ww_mesh = openmc.RegularMesh()
            ww_mesh.dimension = (6, 6, 6)
            ww_mesh.lower_left = (0.0, 0.0, 0.0)
            ww_mesh.upper_right = (30.0, 30.0, 30.0)
            model.settings.weight_window_generators = \
                openmc.WeightWindowGenerator(
                    method="fw_cadis", mesh=ww_mesh, max_realizations=40)
        model.settings.particles = 10
        model.settings.inactive = 20
        model.settings.batches = 40

        harness = MGXSTestHarness('statepoint.40.h5', model)
        harness.main()
