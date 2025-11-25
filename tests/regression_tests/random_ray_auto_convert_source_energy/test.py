import os

import openmc
from openmc.examples import pwr_pin_cell
from openmc import RegularMesh
from openmc.utility_funcs import change_directory
import pytest

from tests.testing_harness import TolerantPyAPITestHarness


class MGXSTestHarness(TolerantPyAPITestHarness):
    def _cleanup(self):
        super()._cleanup()
        f = 'mgxs.h5'
        if os.path.exists(f):
            os.remove(f)


@pytest.mark.parametrize("method", ["stochastic_slab", "infinite_medium"])
def test_random_ray_auto_convert_source_energy(method):
    with change_directory(method):
        openmc.reset_auto_ids()

        # Start with a normal continuous energy model
        model = pwr_pin_cell()

        # Define the source energy distribution. Since this is a PWR pin
        # cell, we use the default Watt fission spectrum.
        source_energy = openmc.stats.Watt()

        # Convert to a multi-group model
        model.convert_to_multigroup(
            method=method, groups='CASMO-8', nparticles=100,
            overwrite_mgxs_library=False, mgxs_path="mgxs.h5",
            source_energy=source_energy
        )

        # Convert to a random ray model
        model.convert_to_random_ray()

        # Set the number of particles
        model.settings.particles = 100

        # Overlay a basic 2x2 mesh
        n = 2
        mesh = RegularMesh()
        mesh.dimension = (n, n)
        bbox = model.geometry.bounding_box
        mesh.lower_left = (bbox.lower_left[0], bbox.lower_left[1])
        mesh.upper_right = (bbox.upper_right[0], bbox.upper_right[1])
        model.settings.random_ray['source_region_meshes'] = [
            (mesh, [model.geometry.root_universe])]

        # Set the source shape to linear
        model.settings.random_ray['source_shape'] = 'linear'

        harness = MGXSTestHarness('statepoint.10.h5', model)
        harness.main()
