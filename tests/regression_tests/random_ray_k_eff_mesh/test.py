import os

import openmc
from openmc.examples import random_ray_lattice

from tests.testing_harness import TolerantPyAPITestHarness


class MGXSTestHarness(TolerantPyAPITestHarness):
    def _cleanup(self):
        super()._cleanup()
        f = 'mgxs.h5'
        if os.path.exists(f):
            os.remove(f)


def test_random_ray_basic():
    model = random_ray_lattice()

    # The model already has some geometrical subdivisions
    # up to a 10x10 grid in the moderator region. So, we
    # double thee re
    pitch = 1.26
    dim = 40
    mesh = openmc.RegularMesh()
    mesh.dimension = (dim, dim)
    mesh.lower_left = (-pitch, -pitch)
    mesh.upper_right = (pitch, pitch)
    
    root = model.geometry.root_universe
    
    model.settings.random_ray['source_region_meshes'] = [(mesh, [root])]

    harness = MGXSTestHarness('statepoint.10.h5', model)
    harness.main()
