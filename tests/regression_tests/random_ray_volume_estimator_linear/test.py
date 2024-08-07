import os

import openmc
from openmc.utility_funcs import change_directory
from openmc.examples import random_ray_three_region_cube
import pytest

from tests.testing_harness import PyAPITestHarness


class MGXSTestHarness(PyAPITestHarness):
    def _cleanup(self):
        super()._cleanup()
        f = 'mgxs.h5'
        if os.path.exists(f):
            os.remove(f)


@pytest.mark.parametrize("estimator", ["hybrid",
                                       "simulation_averaged",
                                       "naive"
                                       ])
def test_random_ray_volume_estimator_linear(estimator):
    with change_directory(estimator):
        openmc.reset_auto_ids()
        model = random_ray_three_region_cube()
        model.settings.random_ray['source_shape'] = 'linear'
        model.settings.random_ray['volume_estimator'] = estimator
        model.settings.inactive = 20
        model.settings.batches = 40
        harness = MGXSTestHarness('statepoint.40.h5', model)
        harness.main()
