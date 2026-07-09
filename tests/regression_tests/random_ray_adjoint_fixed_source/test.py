import os

import openmc
from openmc.utility_funcs import change_directory
from openmc.examples import random_ray_three_region_cube

from tests.testing_harness import TolerantPyAPITestHarness


class MGXSTestHarness(TolerantPyAPITestHarness):
    def _cleanup(self):
        super()._cleanup()
        f = 'mgxs.h5'
        if os.path.exists(f):
            os.remove(f)


def test_random_ray_adjoint_fixed_source():
    model = random_ray_three_region_cube()
    model.settings.random_ray['adjoint'] = True
    model.settings.random_ray['volume_estimator'] = 'naive'
    model.settings.particles = 500
    harness = MGXSTestHarness('statepoint.10.h5', model)
    harness.main()


def test_random_ray_adjoint_fixed_source_adaptive_starved():
    # Ray-starved adaptive case (~20% forward / ~40% adjoint miss rate): the
    # adjoint (second) solve makes real end-of-inactive demotion decisions on
    # its own accumulated flux, guarding both the adaptive machinery in
    # adjoint mode and the clearing of the forward solve's accumulated flux
    # between the two solves (which would otherwise swamp the demotion test).
    with change_directory('adaptive_starved'):
        openmc.reset_auto_ids()
        model = random_ray_three_region_cube()
        model.settings.random_ray['adjoint'] = True
        model.settings.random_ray['volume_estimator'] = 'adaptive'
        model.settings.particles = 10
        model.settings.inactive = 20
        model.settings.batches = 40
        harness = MGXSTestHarness('statepoint.40.h5', model)
        harness.main()
