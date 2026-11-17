"""Runtime spectrum access uses native NNDC U235 prompt-fission data."""
from copy import deepcopy

import h5py
import numpy as np
import pytest

import openmc
import openmc.lib
from openmc import exceptions as exc


@pytest.fixture
def loaded_spectrum(sphere_model, run_in_tmpdir, mpi_intracomm):
    try:
        sphere_model.init_lib(output=False, intracomm=mpi_intracomm)
        index = openmc.lib.nuclides['U235']._index
        yield index, openmc.lib.get_emission_spectrum(index, 18)
    finally:
        sphere_model.finalize_lib()


@pytest.fixture
def initialized_spectrum(loaded_spectrum):
    openmc.lib.simulation_init()
    return loaded_spectrum


def assert_spectrum_equal(actual, expected):
    np.testing.assert_array_equal(actual['incident_energy'],
                                  expected['incident_energy'])
    assert len(actual['tables']) == len(expected['tables'])
    for a, b in zip(actual['tables'], expected['tables']):
        assert a['n_discrete'] == b['n_discrete']
        assert a['interpolation'] == b['interpolation']
        for key in ('energy', 'pdf', 'cdf'):
            np.testing.assert_array_equal(a[key], b[key])


def test_read_native_spectrum(loaded_spectrum):
    index, original = loaded_spectrum
    library = openmc.data.DataLibrary.from_xml()
    path = library.get_by_material('U235')['path']
    with h5py.File(path) as f:
        group = f['U235/reactions/reaction_018/product_0/distribution_0/energy']
        np.testing.assert_array_equal(original['incident_energy'],
                                      group['energy'][()])
        data = group['distribution']
        offsets = list(data.attrs['offsets']) + [data.shape[1]]
        for i, table in enumerate(original['tables']):
            np.testing.assert_array_equal(
                [table['energy'], table['pdf'], table['cdf']],
                data[:, offsets[i]:offsets[i + 1]])
            assert table['n_discrete'] == data.attrs['n_discrete_lines'][i]
            assert table['interpolation'] == data.attrs['interpolation'][i]
    saved = deepcopy(original)
    original['tables'][0]['pdf'][:] = 0.0
    assert_spectrum_equal(openmc.lib.get_emission_spectrum(index, 18), saved)


@pytest.mark.parametrize('interpolation', [1, 2])
def test_replace_and_read_back(initialized_spectrum, interpolation):
    index, original = initialized_spectrum
    modified = deepcopy(original)
    for table in modified['tables']:
        table.update(energy=np.array([0.0, 1.e6, 2.e6]),
                     pdf=np.array([1.0, 3.0, 1.0]),
                     cdf=np.array([99.0]),  # CDF is deliberately ignored.
                     interpolation=interpolation, n_discrete=0)
    openmc.lib.set_emission_spectrum(index, 18, modified)
    actual = openmc.lib.get_emission_spectrum(index, 18)
    np.testing.assert_array_equal(actual['incident_energy'],
                                  original['incident_energy'])
    for table in actual['tables']:
        np.testing.assert_array_equal(table['energy'], [0.0, 1.e6, 2.e6])
        np.testing.assert_allclose(table['pdf'], np.array([1, 3, 1]) / 4.e6)
        np.testing.assert_allclose(table['cdf'],
                                   [0.0, 0.25 if interpolation == 1 else 0.5, 1.0])
        assert table['interpolation'] == interpolation
    assert not np.array_equal(actual['tables'][0]['pdf'],
                              original['tables'][0]['pdf'])


@pytest.mark.parametrize('case', [
    'short_grid', 'long_grid', 'changed_grid', 'nan_grid', 'short_tables',
    'long_tables', 'short_both', 'long_both', 'empty', 'zero_pdf',
    'negative_pdf', 'nan_pdf', 'infinite_energy', 'negative_energy',
    'unordered_energy', 'duplicate_energy', 'one_continuous_point',
    'negative_discrete', 'excess_discrete', 'bad_interpolation',
    'fractional_interpolation', 'wrong_length', 'scalar_energy',
    'incompatible_discrete_only',
])
def test_invalid_replacement_is_atomic(initialized_spectrum, case):
    index, original = initialized_spectrum
    bad = deepcopy(original)
    # Make the first table valid but different, so a failure in the final table
    # detects partial native mutation rather than merely an exception.
    first = bad['tables'][0]
    first['pdf'] *= np.exp(-first['energy'] / 3.e5)
    last = bad['tables'][-1]
    if case in ('short_grid', 'short_both'):
        bad['incident_energy'] = bad['incident_energy'][:-1]
    if case in ('long_grid', 'long_both'):
        bad['incident_energy'] = np.append(bad['incident_energy'], 3.e7)
    if case == 'changed_grid':
        bad['incident_energy'][0] += 1.0
    if case == 'nan_grid':
        bad['incident_energy'][0] = np.nan
    if case in ('short_tables', 'short_both'):
        bad['tables'].pop()
    if case in ('long_tables', 'long_both'):
        bad['tables'].append(deepcopy(last))
    if case == 'empty':
        last['energy'] = last['pdf'] = []
    if case == 'zero_pdf':
        last['pdf'][:] = 0.0
    if case == 'negative_pdf':
        last['pdf'][0] = -1.0
    if case == 'nan_pdf':
        last['pdf'][0] = np.nan
    if case == 'infinite_energy':
        last['energy'][-1] = np.inf
    if case == 'negative_energy':
        last['energy'][0] = -1.0
    if case == 'unordered_energy':
        last['energy'] = last['energy'][::-1]
    if case == 'duplicate_energy':
        last['energy'][1] = last['energy'][0]
    if case == 'one_continuous_point':
        last['n_discrete'] = len(last['energy']) - 1
    if case == 'negative_discrete':
        last['n_discrete'] = -1
    if case == 'excess_discrete':
        last['n_discrete'] = len(last['energy']) + 1
    if case == 'bad_interpolation':
        last['interpolation'] = 3
    if case == 'fractional_interpolation':
        last['interpolation'] = 1.5
    if case == 'wrong_length':
        last['pdf'] = last['pdf'][:-1]
    if case == 'scalar_energy':
        last['energy'] = 1.0
    if case == 'incompatible_discrete_only':
        last['n_discrete'] = len(last['energy'])
    with pytest.raises((ValueError, exc.DataError)):
        openmc.lib.set_emission_spectrum(index, 18, bad)
    assert_spectrum_equal(openmc.lib.get_emission_spectrum(index, 18), original)


def test_requires_simulation_initialization(loaded_spectrum):
    index, spectrum = loaded_spectrum
    with pytest.raises(exc.AllocationError, match='initialized'):
        openmc.lib.set_emission_spectrum(index, 18, spectrum)
    openmc.lib.simulation_init()
    openmc.lib.set_emission_spectrum(index, 18, spectrum)
    saved = openmc.lib.get_emission_spectrum(index, 18)
    openmc.lib.simulation_finalize()
    with pytest.raises(exc.AllocationError, match='initialized'):
        openmc.lib.set_emission_spectrum(index, 18, spectrum)
    assert_spectrum_equal(openmc.lib.get_emission_spectrum(index, 18), saved)


def test_requires_loaded_library():
    spectrum = {'type': 'continuous-tabular', 'incident_energy': [1.0, 2.0],
                'tables': [{'energy': [0.0, 1.0], 'pdf': [1.0, 1.0],
                            'interpolation': 2, 'n_discrete': 0}] * 2}
    with pytest.raises(exc.AllocationError, match='initialized'):
        openmc.lib.set_emission_spectrum(0, 18, spectrum)


@pytest.mark.parametrize('mt,product,distribution', [
    (2, 0, 0), (999, 0, 0), (18, -1, 0), (18, 0, -1), (18, 0, 999),
])
def test_unsupported_selection(initialized_spectrum, mt, product, distribution):
    index, original = initialized_spectrum
    with pytest.raises(exc.DataError):
        openmc.lib.get_emission_spectrum(index, mt, product, distribution)
    with pytest.raises(exc.DataError):
        openmc.lib.set_emission_spectrum(index, mt, original, product, distribution)
    assert_spectrum_equal(openmc.lib.get_emission_spectrum(index, 18), original)


def test_existing_xs_editing(initialized_spectrum):
    index, _ = initialized_spectrum
    energy = openmc.lib.nuclide_energy_grid(index)
    original = openmc.lib.get_nuclide_xs(index, 18)
    threshold = openmc.lib.get_reaction_threshold_energy(index, 18)
    values = np.zeros_like(energy)
    values[energy >= threshold] = original * 1.01
    openmc.lib.new_nuclide_xs_with_threshold(index, 18, energy, values, threshold)
    openmc.lib.nuclide_rebuild_derived_xs(index)
    np.testing.assert_allclose(openmc.lib.get_nuclide_xs(index, 18),
                               original * 1.01)


@pytest.mark.parametrize('case', ['offset_start', 'offset_end', 'offset_order',
                                 'native_grid', 'native_count'])
def test_native_validation_is_atomic(initialized_spectrum, case):
    index, original = initialized_spectrum
    grid = original['incident_energy'].copy()
    n = len(grid)
    offsets = np.arange(0, 2 * n + 1, 2, dtype=np.int32)
    interpolation = np.full(n, 2, dtype=np.int32)
    discrete = np.zeros(n, dtype=np.int32)
    energy = np.tile([0.0, 2.e6], n)
    pdf = np.tile([1.0, 1.0], n)
    if case == 'offset_start':
        offsets[0] = 1
    elif case == 'offset_end':
        offsets[-1] -= 1
    elif case == 'offset_order':
        offsets[-2] = offsets[-1] + 1
    elif case == 'native_grid':
        grid[-1] += 1.0
    elif case == 'native_count':
        n -= 1
    # Bypass Python validation to check that the C entry point is safe itself.
    with pytest.raises(exc.DataError):
        openmc.lib._dll.openmc_nuclide_set_emission_spectrum(
            index, 18, 0, 0, grid, n, offsets, interpolation, discrete,
            energy, pdf, len(energy))
    assert_spectrum_equal(openmc.lib.get_emission_spectrum(index, 18), original)
