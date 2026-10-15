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


def test_random_ray_source_gradient_limiter():
    # A linear source run with the gradient limiter enabled and firing in
    # both of its regimes: the naive volume estimator and an overlay
    # source-region mesh leave the example's optically thin interior with
    # noisy fitted gradients, and the absorber's steep attenuation over
    # regions a few mean free paths thick gives physically steep ones, which
    # the limiter clips as well. The example's three cubic regions are
    # replaced by spherical ones so that the curved boundaries cut the mesh
    # cells into pieces whose centroids sit off-center in their bounding
    # boxes, which is where the limiter's bound differs from a symmetric one.
    openmc.reset_auto_ids()
    model = random_ray_three_region_cube()
    source_mat, void_mat, absorber_mat = model.materials
    width = 30.0
    x0 = openmc.XPlane(0.0, boundary_type='reflective')
    y0 = openmc.YPlane(0.0, boundary_type='reflective')
    z0 = openmc.ZPlane(0.0, boundary_type='reflective')
    x1 = openmc.XPlane(width, boundary_type='vacuum')
    y1 = openmc.YPlane(width, boundary_type='vacuum')
    z1 = openmc.ZPlane(width, boundary_type='vacuum')
    domain = +x0 & -x1 & +y0 & -y1 & +z0 & -z1
    source_sphere = openmc.Sphere(r=5.0)
    void_sphere = openmc.Sphere(r=12.5)
    model.geometry = openmc.Geometry([
        openmc.Cell(fill=source_mat, region=-source_sphere & domain),
        openmc.Cell(fill=void_mat,
                    region=+source_sphere & -void_sphere & domain),
        openmc.Cell(fill=absorber_mat, region=+void_sphere & domain),
    ])
    model.settings.source[0].constraints = {'domains': [source_mat]}
    model.settings.random_ray['source_shape'] = 'linear'
    model.settings.random_ray['source_gradient_limiter'] = True
    model.settings.random_ray['volume_estimator'] = 'naive'
    mesh = openmc.RegularMesh()
    mesh.lower_left = (0.0, 0.0, 0.0)
    mesh.upper_right = (width, width, width)
    mesh.dimension = (12, 12, 12)
    model.settings.random_ray['source_region_meshes'] = [
        (mesh, [model.geometry.root_universe])]
    model.settings.inactive = 30
    model.settings.batches = 60
    harness = MGXSTestHarness('statepoint.60.h5', model)
    harness.main()
