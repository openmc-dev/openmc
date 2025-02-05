import os

import openmc
from openmc.examples import random_ray_three_region_cube
from openmc.utility_funcs import change_directory
import pytest

from tests.testing_harness import TolerantPyAPITestHarness


class MGXSTestHarness(TolerantPyAPITestHarness):
    def _cleanup(self):
        super()._cleanup()
        f = 'mgxs.h5'
        if os.path.exists(f):
            os.remove(f)

def make_mesh(dim):
    width = 30.0
    mesh = openmc.RegularMesh()
    mesh.dimension = (dim, dim, dim)
    mesh.lower_left = (0.0, 0.0, 0.0)
    mesh.upper_right = (width, width, width)
    return mesh

@pytest.mark.parametrize("shape", ["flat", "linear"])
def test_random_ray_fixed_source_mesh(shape):
    with change_directory(shape):
        openmc.reset_auto_ids()
        model = random_ray_three_region_cube()

        # We will apply three different mesh resolutions to three different domain types
        source_universe = model.geometry.get_universes_by_name('source universe')[0]
        void_cell = model.geometry.get_cells_by_name('infinite void region')[0]
        absorber_mat = model.geometry.get_materials_by_name('absorber')[0]
        
        mesh1 = make_mesh(24)
        mesh2 = make_mesh(36)
        mesh3 = make_mesh(30)

        model.settings.random_ray['source_region_meshes'] = [(mesh1, [source_universe]), (mesh2, [void_cell]), (mesh3, [absorber_mat])]

        # We also test flat/linear source shapes to ensure they are both
        # working correctly with the mesh overlay logic
        model.settings.random_ray['source_shape'] = shape

        harness = MGXSTestHarness('statepoint.10.h5', model)
        harness.main()
