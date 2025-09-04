import os

import openmc
from openmc.utility_funcs import change_directory
from openmc.examples import random_ray_three_region_sphere
import pytest

from tests.testing_harness import TolerantPyAPITestHarness


class MGXSTestHarness(TolerantPyAPITestHarness):
    def _cleanup(self):
        super()._cleanup()
        f = "mgxs.h5"
        if os.path.exists(f):
            os.remove(f)


def test_random_ray_fixed_source():
    openmc.reset_auto_ids()
    model = random_ray_three_region_sphere()
    source = model.settings.source[0]

    dists = []
    probs = []

    dists.append(openmc.stats.Discrete([1e2, 1e6], [0.5, 0.5]))
    probs.append(1.0)

    dists.append(openmc.stats.Watt())
    probs.append(2.0)

    dists.append(openmc.stats.Maxwell(1e6))
    probs.append(1.0)

    dists.append(openmc.stats.Uniform(1, 1e6))
    probs.append(1.0)

    dists.append(openmc.stats.PowerLaw(1, 1e6, -1.0))
    probs.append(1.0)

    source.energy = openmc.stats.Mixture(probs, dists)

    harness = MGXSTestHarness("statepoint.10.h5", model)
    harness.main()
