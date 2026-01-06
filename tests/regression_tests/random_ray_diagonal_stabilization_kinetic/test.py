import os

import openmc
from openmc.examples import pwr_pin_cell
from openmc import RegularMesh
from openmc.utility_funcs import change_directory
import numpy as np
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
def test_random_ray_diagonal_stabilization(time_method):
    with change_directory(time_method):
        openmc.reset_auto_ids()

        # Start with a normal continuous energy model
        model = pwr_pin_cell()

        # Convert to a multi-group model, with 70 group XS
        # and transport correction enabled. This will generate
        # MGXS data with some negatives on the diagonal, in order
        # to trigger diagonal correction.
        model.convert_to_multigroup(
            method='material_wise', energy_groups='CASMO-70', nparticles=30,
            overwrite_mgxs_library=True, mgxs_path="mgxs.h5", correction='P0',
            kinetic=True, num_delayed_groups=6
        )

        model.settings.timestep_parameters['n_timesteps'] = 5
        density_timeseries = np.linspace(1, 0.95, 100)
        model.materials[2].set_density('macro', density=1.0, density_timeseries=density_timeseries)

        # Convert to a random ray model
        model.convert_to_random_ray()

        # Set the number of particles
        model.settings.particles = 100

        # Overlay an 8x8 mesh
        n = 8
        mesh = RegularMesh()
        mesh.dimension = (n, n)
        bbox = model.geometry.bounding_box
        mesh.lower_left = (bbox.lower_left[0], bbox.lower_left[1])
        mesh.upper_right = (bbox.upper_right[0], bbox.upper_right[1])
        model.settings.random_ray['source_region_meshes'] = [
            (mesh, [model.geometry.root_universe])]
        model.settings.random_ray['time_derivative_method'] = time_method

        # Explicitly set the diagonal stabilization rho (default is otherwise 1.0).
        # Note that if we set this to 0.0 (thus distabling stabilization), the
        # problem should fail due to instability, so this is actually a good test
        # problem.
        model.settings.random_ray['diagonal_stabilization_rho'] = 0.5

        # If rho was 0.0, the instability would cause failure after iteration 14,
        # so we go a little past that.
        model.settings.inactive = 15
        model.settings.batches = 20

        harness = KineticMGXSTestHarness('statepoint.20', 6, model)
        harness.main()
