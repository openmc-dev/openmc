"""Angular snapshots validated against native transport sampling."""
from copy import deepcopy

import numpy as np
import pytest
from scipy.integrate import quad

import openmc
import openmc.lib
from openmc import exceptions
from sandbox.examples.runtime_angular_distribution import (
    angular_mass, incident_weights, kalbach_pdf, table_mass, validate_samples)
from tests.unit_tests.test_lib_secondary import NUCLIDES, native_data


@pytest.mark.parametrize('nuclide', NUCLIDES)
def test_angular_normalization(native_data, nuclide):
    _, names = native_data
    if nuclide not in names:
        pytest.skip(f'{nuclide} not installed')
    for product in openmc.lib.get_reaction_products(nuclide):
        for descriptor in product['distributions']:
            law = openmc.lib.get_emission_distribution(
                nuclide, product['mt'], product['product_index'],
                descriptor['distribution_index'])
            if law['type'] == 'uncorrelated':
                angle = law['angle']
                for table in angle.get('tables', []):
                    assert -1 <= table['x'][0] < table['x'][-1] <= 1
                    assert table_mass(table) == pytest.approx(1, abs=1.e-6)
                    assert table['cdf'][0] == 0
                    assert table['cdf'][-1] == 1
            elif law['type'] == 'correlated':
                for table in law['tables']:
                    for angle in table['angle']:
                        assert table_mass(angle) == pytest.approx(1, abs=1.e-6)
            elif law['type'] == 'kalbach-mann':
                for i in {0, len(law['tables']) // 2, len(law['tables']) - 1}:
                    table = law['tables'][i]
                    for j in {0, len(table['a']) // 2, len(table['a']) - 1}:
                        r, a = table['r'][j], table['a'][j]
                        integral = quad(
                            lambda mu: kalbach_pdf(mu, r, a), -1, 1)[0]
                        assert integral == pytest.approx(1)


@pytest.mark.parametrize('mt', [2, 24])
def test_native_angular_sampler(native_data, mt):
    angle = openmc.lib.get_angular_distribution('Li6', mt)
    assert angle['type'] == 'tabular-angle'
    _, _, probability = validate_samples('Li6', mt, 14.1e6, angle)
    assert probability > 1.e-5
    # Between incident points, test the stochastic mixture rather than choosing
    # the nearest grid point. Check reproducibility without modifying RNG streams.
    first = openmc.lib.sample_angular_distribution('Li6', mt, 14.1e6, 100, seed=88)
    second = openmc.lib.sample_angular_distribution('Li6', mt, 14.1e6, 100, seed=88)
    np.testing.assert_array_equal(first, second)
    assert angular_mass(angle, 14.1e6) == pytest.approx(1, abs=1.e-6)


def test_angular_copy(native_data):
    angle = openmc.lib.get_angular_distribution('Li6', 2)
    saved = deepcopy(angle)
    angle['incident_energy'][:] = 0
    for table in angle['tables']:
        for key in ('x', 'pdf', 'cdf'):
            table[key][:] = -123
        table['interpolation'] = -1
    reread = openmc.lib.get_emission_distribution('Li6', 2)['angle']
    np.testing.assert_array_equal(reread['incident_energy'], saved['incident_energy'])
    for actual, expected in zip(reread['tables'], saved['tables']):
        for key in ('x', 'pdf', 'cdf'):
            np.testing.assert_array_equal(actual[key], expected[key])
        assert actual['interpolation'] == expected['interpolation']


@pytest.mark.parametrize('nuclide,mt', [('Be9', 16), ('C12', 5)])
def test_conditional_rejected(native_data, nuclide, mt):
    _, names = native_data
    if nuclide not in names:
        pytest.skip(f'{nuclide} not installed')
    with pytest.raises(ValueError, match='conditional.*get_emission_distribution'):
        openmc.lib.get_angular_distribution(nuclide, mt)
    with pytest.raises(exceptions.InvalidArgumentError, match='conditional'):
        openmc.lib.sample_angular_distribution(nuclide, mt, 14.1e6, 10)


def test_nbody_isotropic(native_data):
    angle = openmc.lib.get_angular_distribution('H2', 16)
    assert angle == {'type': 'isotropic'}
    _, _, probability = validate_samples('H2', 16, 14.1e6, angle)
    assert probability > 1.e-5


@pytest.mark.parametrize('energy', [-1, np.nan, np.inf])
def test_bad_energy(native_data, energy):
    with pytest.raises(exceptions.InvalidArgumentError):
        openmc.lib.sample_angular_distribution('Li6', 2, energy, 10)


@pytest.mark.parametrize('kwargs', [{'n_samples': -1}, {'n_samples': True},
                                   {'seed': -1}, {'seed': 2**64}, {'seed': True}])
def test_bad_sampling_arguments(native_data, kwargs):
    args = dict(n_samples=10)
    args.update(kwargs)
    with pytest.raises((ValueError, TypeError)):
        openmc.lib.sample_angular_distribution('Li6', 2, 14.1e6, **args)


def test_incident_interpolation_and_histogram_integral():
    # This asymmetric histogram would fail a blind trapezoidal integration.
    table = {'x': np.array([-1., 0., 1.]), 'pdf': np.array([.25, .75, 9.]),
             'interpolation': 1}
    assert table_mass(table) == 1
    assert table_mass(table, -.5, .5) == .5
    angle = {'incident_energy': np.array([1., 3.]), 'tables': [table, table]}
    assert incident_weights(angle, 0.) == [(0, 1.)]
    assert incident_weights(angle, 2.) == [(0, .5), (1, .5)]
    assert incident_weights(angle, 4.) == [(1, 1.)]
