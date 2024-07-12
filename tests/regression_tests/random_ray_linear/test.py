import os

import openmc
from openmc.examples import random_ray_lattice
from openmc.utility_funcs import change_directory
import pytest

from tests.testing_harness import TolerantPyAPITestHarness


class MGXSTestHarness(TolerantPyAPITestHarness):
    def _cleanup(self):
        super()._cleanup()
        f = 'mgxs.h5'
        if os.path.exists(f):
            os.remove(f)


@pytest.mark.parametrize("shape", ["linear", "linear_xy"])
def test_random_ray_source(shape):
    with change_directory(shape):
        openmc.reset_auto_ids()
        model = random_ray_lattice()
        model.settings.random_ray['source_shape'] = shape
        harness = MGXSTestHarness('statepoint.10.h5', model)
        harness.main()
