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


# A deliberately ray-starved configuration (~20% source region miss rate).
# The volume estimators only differ meaningfully when regions are missed or
# sparsely hit, so a starved run exercises every estimator code path. At
# this density the per-estimator volume choices, the miss treatments, and
# for the adaptive estimator the strong-source (kappa) demotion, the
# hit-starved demotion, the end-of-inactive converged-negative demotion,
# and the previous-flux miss treatment all fire.
@pytest.mark.parametrize("estimator", ["hybrid",
                                       "simulation_averaged",
                                       "naive",
                                       "adaptive",
                                       "strict_adaptive"
                                       ])
def test_random_ray_volume_estimator(estimator):
    with change_directory(estimator):
        openmc.reset_auto_ids()
        model = random_ray_three_region_cube()
        model.settings.random_ray['volume_estimator'] = estimator
        model.settings.particles = 10
        model.settings.inactive = 20
        model.settings.batches = 40

        harness = MGXSTestHarness('statepoint.40.h5', model)
        harness.main()
