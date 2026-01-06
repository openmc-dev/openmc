import os

import openmc
from openmc.utility_funcs import change_directory
from openmc.examples import random_ray_lattice
import pytest

from tests.testing_harness import KineticTolerantPyAPITestHarness


class KineticMGXSTestHarness(KineticTolerantPyAPITestHarness):
    def _cleanup(self):
        super()._cleanup()
        f = 'mgxs.h5'
        if os.path.exists(f):
            os.remove(f)


@pytest.mark.parametrize("time_method", ["isotropic",
                                         "propagation"])
def test_random_ray_k_eff_mesh(time_method):
    with change_directory(time_method):
        model = random_ray_lattice(kinetic=True)
        model.settings.timestep_parameters['n_timesteps'] = 5
        model.settings.random_ray['time_method'] = time_method
        model.settings.batches = 400
        model.settings.inactive = 200

        # The model already has some geometrical subdivisions
        # up to a 10x10 grid in the moderator region. So, we
        # increase the resolution 40x40 applied over the full
        # 2x2 lattice.
        pitch = 1.26
        dim = 40
        mesh = openmc.RegularMesh()
        mesh.dimension = (dim, dim)
        mesh.lower_left = (-pitch, -pitch)
        mesh.upper_right = (pitch, pitch)

        root = model.geometry.root_universe

        model.settings.random_ray['source_region_meshes'] = [(mesh, [root])]

        harness = KineticMGXSTestHarness('statepoint.400', 6, model)
        harness.main()
