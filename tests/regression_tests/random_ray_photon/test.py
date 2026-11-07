from pathlib import Path

import numpy as np
import openmc
import pytest

from openmc.utility_funcs import change_directory
from tests.regression_tests import config
from tests.testing_harness import TolerantPyAPITestHarness


def _make_model():
    group_edges = [1.0e3, 1.0e5, 2.0e6]
    groups = openmc.mgxs.EnergyGroups(group_edges=group_edges)
    photon = openmc.XSdata('photon', groups)
    photon.order = 0
    photon.set_total([0.2, 0.3])
    photon.set_absorption([0.04, 0.1])
    photon.set_scatter_matrix(np.array([
        [[0.08], [0.06]],
        [[0.0], [0.1]],
    ]))

    mgxs = openmc.MGXSLibrary(groups, particle_type='photon')
    mgxs.add_xsdata(photon)
    mgxs.export_to_hdf5('mgxs.h5')

    material = openmc.Material()
    material.set_density('macro', 1.0)
    material.add_macroscopic('photon')

    boundary = openmc.model.RectangularParallelepiped(
        -5.0, 5.0, -5.0, 5.0, -5.0, 5.0, boundary_type='vacuum')
    cell = openmc.Cell(fill=material, region=-boundary)

    model = openmc.Model()
    model.geometry = openmc.Geometry([cell])
    model.materials = openmc.Materials([material])
    model.materials.cross_sections = 'mgxs.h5'

    mesh = openmc.RegularMesh()
    mesh.dimension = (2, 1, 1)
    mesh.lower_left = (-5.0, -5.0, -5.0)
    mesh.upper_right = (5.0, 5.0, 5.0)

    tally = openmc.Tally()
    tally.filters = [
        openmc.MeshFilter(mesh),
        openmc.EnergyFilter(group_edges),
        openmc.ParticleFilter('photon'),
    ]
    tally.scores = ['flux', 'total']
    model.tallies.append(tally)

    model.settings.energy_mode = 'multi-group'
    model.settings.run_mode = 'fixed source'
    model.settings.photon_transport = True
    model.settings.batches = 10
    model.settings.inactive = 3
    model.settings.particles = 200
    model.settings.source = openmc.IndependentSource(
        particle='photon',
        space=openmc.stats.Point((-2.5, 0.0, 0.0)),
        energy=openmc.stats.Discrete([1.0e6], [1.0]))
    model.settings.random_ray = {
        'distance_inactive': 50.0,
        'distance_active': 200.0,
        'ray_source': openmc.IndependentSource(
            particle='photon',
            space=openmc.stats.Box(
                (-5.0, -5.0, -5.0), (5.0, 5.0, 5.0))),
        'source_region_meshes': [(mesh, [model.geometry.root_universe])],
    }

    return model


@pytest.fixture
def model():
    return _make_model()


class PhotonRandomRayTestHarness(TolerantPyAPITestHarness):
    def _get_results(self, hash_output=False):
        with openmc.StatePoint(self._sp_name) as statepoint:
            tally = statepoint.tallies[1]
            if not np.all(tally.mean > 0.0):
                raise AssertionError(
                    'Photon random ray tally scores must be positive')
        return super()._get_results(hash_output)

    def _cleanup(self):
        super()._cleanup()
        Path('mgxs.h5').unlink(missing_ok=True)


def test_random_ray_photon(model):
    harness = PhotonRandomRayTestHarness('statepoint.10.h5', model)
    harness.main()


@pytest.mark.parametrize(
    ('invalid_input', 'error'),
    [
        ('source mismatch', 'does not match'),
        ('eigenvalue', 'only supported in fixed-source mode'),
        ('adjoint', 'Adjoint photon'),
    ],
)
def test_random_ray_photon_validation(tmp_path, invalid_input, error):
    with change_directory(tmp_path):
        model = _make_model()
        if invalid_input == 'source mismatch':
            model.settings.source = openmc.IndependentSource(
                particle='neutron',
                space=openmc.stats.Point((-2.5, 0.0, 0.0)),
                energy=openmc.stats.Discrete([1.0e6], [1.0]))
        elif invalid_input == 'eigenvalue':
            model.settings.run_mode = 'eigenvalue'
        else:
            model.settings.random_ray['adjoint'] = True

        with pytest.raises(RuntimeError, match=error):
            model.run(openmc_exec=config['exe'])
