import os

import openmc
from openmc.examples import pwr_pin_cell
from openmc import RegularMesh
from openmc.utility_funcs import change_directory
import pytest
import numpy as np

from tests.testing_harness import KineticTolerantPyAPITestHarness


class KineticMGXSTestHarness(KineticTolerantPyAPITestHarness):
    def _cleanup(self):
        super()._cleanup()
        f = 'mgxs.h5'
        if os.path.exists(f):
            os.remove(f)


@pytest.mark.parametrize("generation_method, time_method", [("material_wise", "isotropic"),
                                                            ("stochastic_slab",
                                                             "isotropic"),
                                                            ("infinite_medium",
                                                             "isotropic"),
                                                            ("material_wise",
                                                             "propagation"),
                                                            ("stochastic_slab",
                                                             "propagation"),
                                                            ("infinite_medium",
                                                             "propagation"),
                                                            ])
def test_random_ray_auto_convert(generation_method, time_method):
    with change_directory(f'{generation_method}/{time_method}'):
        openmc.reset_auto_ids()

        # Start with a normal continuous energy model
        model = pwr_pin_cell()

        # Convert to a multi-group model
        model.convert_to_multigroup(
            method=generation_method, energy_groups='CASMO-2', nparticles=30,
            overwrite_mgxs_library=False, mgxs_path="mgxs.h5", kinetic=True,
            num_delayed_groups=6
        )

        model.settings.timestep_parameters['n_timesteps'] = 5
        density_timeseries = np.linspace(1, 0.95, 100)
        model.materials[2].set_density(
            'macro', density=1.0, density_timeseries=density_timeseries)

        # Convert to a random ray model
        model.convert_to_random_ray()

        # Set the number of particles
        model.settings.particles = 100

        # Overlay an basic 8x8 mesh
        n = 8
        mesh = RegularMesh()
        mesh.dimension = (n, n)
        bbox = model.geometry.bounding_box
        mesh.lower_left = (bbox.lower_left[0], bbox.lower_left[1])
        mesh.upper_right = (bbox.upper_right[0], bbox.upper_right[1])
        model.settings.random_ray['source_region_meshes'] = [
            (mesh, [model.geometry.root_universe])]
        model.settings.random_ray['time_derivative_method'] = time_method

        harness = KineticMGXSTestHarness("statepoint.10", 6, model)
        harness.main()
