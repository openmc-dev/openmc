import os

import numpy as np
import openmc
from openmc.examples import random_ray_lattice

from tests.testing_harness import TolerantPyAPITestHarness

class MGXSTestHarness(TolerantPyAPITestHarness):
    def _cleanup(self):
        super()._cleanup()
        f = 'mgxs.h5'
        if os.path.exists(f):
            os.remove(f)

def test_random_ray_vacuum():
    model = random_ray_lattice()

    # Convert reflective surfaces to vacuum   
    surfaces = model.geometry.get_all_surfaces()
    for key, surface in surfaces.items():
        if surface.boundary_type == 'reflective':
            surface.boundary_type = 'vacuum'

    # Create a geometry with the two cells and export to XML
    harness = MGXSTestHarness('statepoint.10.h5', model)
    harness.main()