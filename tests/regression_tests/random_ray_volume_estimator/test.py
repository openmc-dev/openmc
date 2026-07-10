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


@pytest.mark.parametrize("estimator", ["hybrid",
                                       "simulation_averaged",
                                       "naive",
                                       "adaptive"
                                       ])
def test_random_ray_volume_estimator(estimator):
    with change_directory(estimator):
        openmc.reset_auto_ids()
        model = random_ray_three_region_cube()
        model.settings.random_ray['volume_estimator'] = estimator

        harness = MGXSTestHarness('statepoint.10.h5', model)
        harness.main()


def test_random_ray_volume_estimator_starved():
    # Deliberately ray-starved adaptive configuration (~20% miss rate): unlike
    # the nominal cases above, this exercises every adaptive mechanism at once
    # -- the strong-source (kappa) demotion, the hit-starved demotion, the
    # end-of-inactive converged-negative demotion, and the previous-flux miss
    # treatment -- so changes to any of those code paths show up in its
    # reference results.
    with change_directory('adaptive_starved'):
        openmc.reset_auto_ids()
        model = random_ray_three_region_cube()
        model.settings.random_ray['volume_estimator'] = 'adaptive'
        model.settings.particles = 10
        model.settings.inactive = 20
        model.settings.batches = 40

        harness = MGXSTestHarness('statepoint.40.h5', model)
        harness.main()
