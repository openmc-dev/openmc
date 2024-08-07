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


@pytest.mark.parametrize("normalize", [True, False])
def test_random_ray_fixed_source(normalize):
    with change_directory(str(normalize)):
        openmc.reset_auto_ids()
        model = random_ray_three_region_cube()
        model.settings.random_ray['volume_normalized_flux_tallies'] = normalize

        harness = MGXSTestHarness('statepoint.10.h5', model)
        harness.main()
