import os
import openmc

from openmc.examples import random_ray_three_region_cube_with_detectors

from tests.testing_harness import TolerantPyAPITestHarness


class MGXSTestHarness(TolerantPyAPITestHarness):
    def _cleanup(self):
        super()._cleanup()
        f = 'mgxs.h5'
        if os.path.exists(f):
            os.remove(f)


def test_random_ray_adjoint_local():
    model = random_ray_three_region_cube_with_detectors()

    detector1_cells = model.geometry.get_cells_by_name("detector 1")
    detector2_cells = model.geometry.get_cells_by_name("detector 2")
    detector_cells = detector1_cells + detector2_cells

    strengths = [1.0]
    midpoints = [100.0]
    energy_distribution = openmc.stats.Discrete(x=midpoints, p=strengths)

    adj_source = openmc.IndependentSource(energy=energy_distribution, constraints={
                                          'domains': detector_cells}, strength=3.14)
    
    model.settings.random_ray['adjoint'] = True
    model.settings.random_ray['adjoint_source'] = adj_source
    model.settings.random_ray['volume_estimator'] = 'naive'
    harness = MGXSTestHarness('statepoint.10.h5', model)
    harness.main()
