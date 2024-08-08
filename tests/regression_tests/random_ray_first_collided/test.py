import os

import numpy as np
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


@pytest.mark.parametrize("shape", ["flat", "linear"])
def test_random_ray_first_collided(shape):
    with change_directory(shape):
        openmc.reset_auto_ids()
        model = random_ray_three_region_cube()
        model.settings.random_ray['source_shape'] = shape
        model.settings.random_ray['first_collided_source'] = True
        strengths = [1.0] # Good - fast group appears largest (besides most thermal)
        midpoints = [100.0]
        energy_dist = openmc.stats.Discrete(x=midpoints,p=strengths)
        lower_left_src = [0.0, 0.0, 0.0]
        upper_right_src = [5.0, 5.0, 5.0]
        spatial_distribution = openmc.stats.Box(lower_left_src, upper_right_src, only_fissionable=False)
        source = openmc.IndependentSource(space=spatial_distribution, energy=energy_dist, strength = 3.14)
        model.settings.source = [source]
        harness = MGXSTestHarness('statepoint.10.h5', model)
        harness.main()