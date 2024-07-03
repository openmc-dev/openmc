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
                                       "segment_corrected"
                                       ])
def test_random_ray_volume_estimator(estimator):
    with change_directory(estimator):
        openmc.reset_auto_ids()
        model = random_ray_three_region_cube()
        model.settings.random_ray['volume_estimator'] = estimator

        harness = MGXSTestHarness('statepoint.10.h5', model)
        harness.main()
