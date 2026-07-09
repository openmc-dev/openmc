import os

import openmc
from openmc.utility_funcs import change_directory
from openmc.examples import random_ray_lattice

from tests.testing_harness import TolerantPyAPITestHarness


class MGXSTestHarness(TolerantPyAPITestHarness):
    def _cleanup(self):
        super()._cleanup()
        f = 'mgxs.h5'
        if os.path.exists(f):
            os.remove(f)


def test_random_ray_k_eff_mesh():
    model = random_ray_lattice()

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

    model.settings.random_ray['volume_estimator'] = 'hybrid'
    harness = MGXSTestHarness('statepoint.10.h5', model)
    harness.main()


def test_random_ray_k_eff_mesh_adaptive_starved():
    # Ray-starved adaptive eigenvalue case (~1.2% miss rate): the subdivision
    # mesh plus a low ray count engages the adaptive demotion machinery
    # (strong-source and hit-starved regions, plus the end-of-inactive
    # demotion decision and the previous-flux miss treatment) in eigenvalue
    # mode, which the nominal 0%-miss eigenvalue tests never exercise.
    with change_directory('adaptive_starved'):
        openmc.reset_auto_ids()
        model = random_ray_lattice()
        pitch = 1.26
        mesh = openmc.RegularMesh()
        mesh.dimension = (40, 40)
        mesh.lower_left = (-pitch, -pitch)
        mesh.upper_right = (pitch, pitch)
        root = model.geometry.root_universe
        model.settings.random_ray['source_region_meshes'] = [(mesh, [root])]
        model.settings.random_ray['volume_estimator'] = 'adaptive'
        model.settings.particles = 25
        model.settings.inactive = 20
        model.settings.batches = 40
        harness = MGXSTestHarness('statepoint.40.h5', model)
        harness.main()
