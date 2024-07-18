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


def test_random_ray_first_collision_rays():
    openmc.reset_auto_ids()
    model = random_ray_three_region_cube()
    model.settings.random_ray['first_collision_source'] = True
    model.settings.random_ray['first_collision_rays'] = 1246
    model.settings.random_ray['first_collision_volume_rays'] = 2000
    strengths = [1.0] 
    midpoints = [100.0]
    energy_dist = openmc.stats.Discrete(x=midpoints,p=strengths)
    lower_left_src = [0.0, 0.0, 0.0]
    upper_right_src = [5.0, 5.0, 5.0]
    spatial_distribution = openmc.stats.Box(lower_left_src, upper_right_src, only_fissionable=False)
    source = openmc.IndependentSource(space=spatial_distribution, energy=energy_dist, strength = 3.14)
    model.settings.source = [source]
    harness = MGXSTestHarness('statepoint.10.h5', model)
    harness.main()