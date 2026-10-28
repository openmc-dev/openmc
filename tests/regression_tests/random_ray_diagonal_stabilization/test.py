import os

import numpy as np
import openmc
import openmc.mgxs
from openmc.examples import pwr_pin_cell
from openmc.utility_funcs import change_directory
from openmc import RegularMesh
import pytest

from tests.testing_harness import TolerantPyAPITestHarness


class MGXSTestHarness(TolerantPyAPITestHarness):
    def _cleanup(self):
        super()._cleanup()
        f = 'mgxs.h5'
        if os.path.exists(f):
            os.remove(f)


def _build_model():
    # Start with a normal continuous energy model
    model = pwr_pin_cell()

    # Convert to a multi-group model, with 70 group XS
    # and transport correction enabled. This will generate
    # MGXS data with some negatives on the diagonal, in order
    # to trigger diagonal correction.
    model.convert_to_multigroup(
        method='material_wise', groups='CASMO-70',
        particles=13,
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

    return model


def _build_homogeneous_model():
    # A homogeneous two-group eigenvalue problem with a negative within-group
    # scattering cross section in the fast group, as a transport correction
    # produces. With chi = (1, 0) and no upscatter the infinite-medium
    # eigenvalue is
    #     k = nu_sigma_f2 * sigma_s12 / ((sigma_t2 - sigma_s22) *
    #                                    (sigma_t1 - sigma_s11))
    #       = 0.375 * 1.0 / (0.5 * 1.5) = 0.5
    # exactly, and the iteration approaches it geometrically from above
    # (1.5, 1.0, 0.75, 0.625, ...). Only a few batches are run, enough for
    # the stored reference to pin that direction.
    groups = openmc.mgxs.EnergyGroups(group_edges=[1e-5, 1.0e3, 20.0e6])
    xs = openmc.XSdata('mat', groups)
    xs.order = 0
    xs.set_total([1.0, 1.0])
    xs.set_absorption([0.5, 0.5])
    xs.set_scatter_matrix(np.array([[[-0.5], [1.0]],
                                    [[0.0], [0.5]]]))
    xs.set_fission([0.0, 0.15])
    xs.set_nu_fission([0.0, 0.375])
    xs.set_chi([1.0, 0.0])
    lib = openmc.MGXSLibrary(groups)
    lib.add_xsdatas([xs])
    lib.export_to_hdf5('mgxs.h5')

    mat = openmc.Material(name='mat')
    mat.set_density('macro', 1.0)
    mat.add_macroscopic(openmc.Macroscopic('mat'))
    model = openmc.Model()
    model.materials = openmc.Materials([mat])
    model.materials.cross_sections = 'mgxs.h5'
    box = openmc.model.RectangularParallelepiped(
        0.0, 10.0, 0.0, 10.0, 0.0, 10.0, boundary_type='reflective')
    cell = openmc.Cell(fill=mat, region=-box)
    model.geometry = openmc.Geometry([cell])

    mesh = RegularMesh()
    mesh.lower_left = (0.0, 0.0, 0.0)
    mesh.upper_right = (10.0, 10.0, 10.0)
    mesh.dimension = (4, 4, 4)

    settings = model.settings
    settings.energy_mode = 'multi-group'
    settings.run_mode = 'eigenvalue'
    settings.particles = 100
    settings.inactive = 3
    settings.batches = 8
    settings.random_ray = {
        'distance_inactive': 30.0,
        'distance_active': 200.0,
        'ray_source': openmc.IndependentSource(
            space=openmc.stats.Box((0.0, 0.0, 0.0), (10.0, 10.0, 10.0))),
        'source_region_meshes': [(mesh, [model.geometry.root_universe])],
    }
    return model


# The transport-corrected (P0) library's negative within-group scattering
# drives some reduced sources negative, which the adaptive estimator must
# handle through its negative-source (strong) treatment and its
# end-of-inactive demotion. The adaptive case pins that interplay.
@pytest.mark.parametrize("estimator", ["hybrid", "adaptive"])
def test_random_ray_diagonal_stabilization(estimator):
    with change_directory(estimator):
        openmc.reset_auto_ids()
        model = _build_model()
        model.settings.random_ray['volume_estimator'] = estimator
        harness = MGXSTestHarness('statepoint.20.h5', model)
        harness.main()


# Every estimator must descend toward the analytic eigenvalue of the
# homogeneous problem. The strict estimator's non-negativity fixup must
# assess the stabilized flux iterate: the negative within-group scattering
# drives the raw fast-group iterate negative early on, which the
# stabilization maps to a positive value. Flooring the raw value first would
# freeze the iteration at the previous iterate, and this problem then climbs
# toward a spurious k of 2.0 instead.
@pytest.mark.parametrize("estimator", ["hybrid", "adaptive", "strict_adaptive"])
def test_random_ray_diagonal_stabilization_homogeneous(estimator):
    with change_directory(f'homogeneous_{estimator}'):
        openmc.reset_auto_ids()
        model = _build_homogeneous_model()
        model.settings.random_ray['volume_estimator'] = estimator
        harness = MGXSTestHarness('statepoint.8.h5', model)
        harness.main()
