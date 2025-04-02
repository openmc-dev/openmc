import os

import openmc
from openmc.examples import random_ray_three_region_cube

from tests.testing_harness import TolerantPyAPITestHarness


class MGXSTestHarness(TolerantPyAPITestHarness):
    def _cleanup(self):
        super()._cleanup()
        f = 'mgxs.h5'
        if os.path.exists(f):
            os.remove(f)


def test_random_ray_point_source_locator():
    model = random_ray_three_region_cube()

    # Overlay subdivided SR mesh to reduce resolution from 2.5cm -> 1cm
    width = 30.0
    mesh = openmc.RegularMesh()
    mesh.dimension = (30, 30, 30)
    mesh.lower_left = (0.0, 0.0, 0.0)
    mesh.upper_right = (width, width, width)
    model.settings.random_ray['source_region_meshes'] = [
        (mesh, [model.geometry.root_universe]),
    ]

    # Define a point source
    strengths = [1.0]
    midpoints = [100.0]
    energy_distribution = openmc.stats.Discrete(x=midpoints, p=strengths)
    spatial_distribution = openmc.stats.Point([2.5, 2.5, 2.5])
    source = openmc.IndependentSource(
        energy=energy_distribution, space=spatial_distribution, strength=3.14)
    model.settings.source = [source]

    # Settings
    model.settings.inactive = 15
    model.settings.batches = 30

    harness = MGXSTestHarness('statepoint.30.h5', model)
    harness.main()
