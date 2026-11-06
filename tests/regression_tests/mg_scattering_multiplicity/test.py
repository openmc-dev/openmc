"""Analytical fixed-source checks for multigroup scattering multiplicity."""

import numpy as np
import pytest

import openmc
from openmc.examples import slab_mg
from tests.regression_tests import config


@pytest.mark.parametrize('scatter_format', ['legendre', 'histogram', 'tabular'])
@pytest.mark.parametrize('multiplicity', [
    pytest.param(np.array([[1.0, 2.0], [1.0, 1.0]]), id='variable'),
    pytest.param(np.ones((2, 2)), id='unit'),
    pytest.param(np.full((2, 2), 1.5), id='constant-1.5'),
])
@pytest.mark.parametrize('survival_biasing', [False, True],
                         ids=['analog', 'survival'])
def test_scattering_multiplicity(run_in_tmpdir, scatter_format, multiplicity,
                                survival_biasing):
    # Matrices use incoming groups as rows, ordered from high to low energy.
    total = np.ones(2)
    absorption = np.array([0.6, 0.8])
    scattering = np.array([[0.2, 0.2], [0.0, 0.2]])
    production = scattering * multiplicity
    groups = openmc.mgxs.EnergyGroups([0.0, 1.0, 2.0])
    data = openmc.XSdata('mat_1', groups)
    data.order = 0
    data.set_total(total)
    data.set_absorption(absorption)
    data.set_scatter_matrix(production[:, :, np.newaxis])
    data.set_multiplicity_matrix(multiplicity)
    if scatter_format != 'legendre':
        data = data.convert_scatter_format(scatter_format, 2)
    library = openmc.MGXSLibrary(groups)
    library.add_xsdata(data)
    library.export_to_hdf5('mgxs.h5')

    # Reuse the multigroup example, replacing its leaking slab by a unit-volume
    # reflected box so that the spatially integrated balance has no leakage.
    model = slab_mg(mgxslib_name='mgxs.h5')
    box = openmc.model.RectangularParallelepiped(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, boundary_type='reflective')
    cell = next(iter(model.geometry.root_universe.cells.values()))
    cell.region = -box
    model.geometry = openmc.Geometry([cell])
    settings = model.settings
    settings.run_mode = 'fixed source'
    settings.inactive = 0
    settings.batches = 20
    settings.particles = 5000
    settings.seed = 12345
    settings.survival_biasing = survival_biasing
    settings.tabular_legendre = {'enable': False}
    settings.source = openmc.IndependentSource(
        space=openmc.stats.Box([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
        energy=openmc.stats.Discrete([1.5], [1.0]))
    settings.output = {'summary': False, 'tallies': False}

    tally = openmc.Tally(name='group balance')
    tally.filters = [openmc.EnergyFilter(groups.group_edges)]
    tally.scores = ['flux', 'absorption']
    tally.estimator = 'tracklength'
    model.tallies = [tally]
    model.export_to_model_xml()
    if config['build_inputs']:
        return

    run_args = {'openmc_exec': config['exe'], 'event_based': config['event']}
    if config['mpi']:
        run_args['mpi_args'] = [config['mpiexec'], '-n', config['mpi_np']]
    openmc.run(**run_args)

    # Unit source in the high-energy group: (diag(total) - B.T) Phi = Q.
    # In the variable-multiplicity case, Phi = (1.25, 0.625), and the
    # integrated absorption rate is 0.6*1.25 + 0.8*0.625 = 1.25.
    expected_flux = np.linalg.solve(np.diag(total) - production.T, [1.0, 0.0])
    if multiplicity[0, 1] == 2.0:
        np.testing.assert_allclose(expected_flux, [1.25, 0.625])
    expected_absorption = absorption * expected_flux
    with openmc.StatePoint('statepoint.20.h5', autolink=False) as statepoint:
        result = statepoint.get_tally(name='group balance')
        for score, expected in [('flux', expected_flux),
                                ('absorption', expected_absorption)]:
            # EnergyFilter bins are in increasing energy; reverse them to the
            # high-to-low convention used by the input cross sections.
            mean = result.get_values(scores=[score]).ravel()[::-1]
            std_dev = result.get_values(scores=[score],
                                        value='std_dev').ravel()[::-1]
            assert np.all(np.isfinite(mean))
            assert np.all(np.isfinite(std_dev))
            # Require enough precision that a large error cannot be hidden by
            # a broad confidence interval. The original variable-multiplicity
            # defect biases the low-energy flux by about 23%, versus at most
            # 9% allowed here. The 0.5% floor accommodates near-zero variance.
            assert np.all(std_dev < 0.015 * expected), (score, std_dev, expected)
            tolerance = np.maximum(6.0 * std_dev, 0.005 * expected)
            assert np.all(np.abs(mean - expected) <= tolerance), (
                score, mean, expected, std_dev)
