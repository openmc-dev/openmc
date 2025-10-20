import os

from openmc.examples import pwr_pin_cell
from openmc import RegularMesh

from tests.testing_harness import TolerantPyAPITestHarness


class MGXSTestHarness(TolerantPyAPITestHarness):
    def _cleanup(self):
        super()._cleanup()
        f = 'mgxs.h5'
        if os.path.exists(f):
            os.remove(f)


def test_random_ray_diagonal_stabilization():
    # Start with a normal continuous energy model
    model = pwr_pin_cell()

    # Convert to a multi-group model, with 70 group XS
    # and transport correction enabled. This will generate
    # MGXS data with some negatives on the diagonal, in order
    # to trigger diagonal correction.
    model.convert_to_multigroup(
        method='material_wise', groups='CASMO-70', nparticles=15,
        overwrite_mgxs_library=True, mgxs_path="mgxs.h5", correction='P0'
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

    # Explicitly set the diagonal stabilization rho (default is otherwise 1.0).
    # Note that if we set this to 0.0 (thus distabling stabilization), the
    # problem should fail due to instability, so this is actually a good test
    # problem.
    model.settings.random_ray['diagonal_stabilization_rho'] = 0.5

    # If rho was 0.0, the instability would cause failure after iteration 14,
    # so we go a little past that.
    model.settings.inactive = 15
    model.settings.batches = 20

    harness = MGXSTestHarness('statepoint.20.h5', model)
    harness.main()
