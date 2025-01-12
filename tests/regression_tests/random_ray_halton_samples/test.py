import os

from openmc.examples import random_ray_lattice

from tests.testing_harness import TolerantPyAPITestHarness


class MGXSTestHarness(TolerantPyAPITestHarness):
    def _cleanup(self):
        super()._cleanup()
        f = 'mgxs.h5'
        if os.path.exists(f):
            os.remove(f)

def test_random_ray_halton_samples():
    model = random_ray_lattice()
    model.settings.random_ray['sample_method'] = 'halton'
    harness = MGXSTestHarness('statepoint.10.h5', model)
    harness.main()
