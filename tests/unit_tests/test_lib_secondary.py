"""Native secondary snapshots compared with the installed evaluated data."""
from copy import deepcopy

import h5py
import numpy as np
import pytest

import openmc
import openmc.lib
from openmc import exceptions


NUCLIDES = ['H1', 'H2', 'He4', 'Li6', 'Li7', 'Be9', 'C12', 'O16',
            'Al27', 'Fe56', 'U235', 'U238', 'Ra223']


def decoded(value):
    return value.decode() if isinstance(value, bytes) else value


def assert_function(actual, dataset):
    kind = decoded(dataset.attrs['type'])
    if kind == 'Polynomial':
        assert actual['type'] == 'polynomial'
        np.testing.assert_array_equal(actual['coefficients'], dataset[()])
    else:
        assert kind == 'Tabulated1D'
        assert actual['type'] == 'tabulated-1d'
        np.testing.assert_array_equal(actual['x'], dataset[0])
        np.testing.assert_array_equal(actual['y'], dataset[1])
        for key in ('breakpoints', 'interpolation'):
            expected = dataset.attrs[key]
            if key == 'interpolation':
                expected = expected % 10
            np.testing.assert_array_equal(actual[key], expected)


def assert_mu(actual, values, interpolation):
    assert actual['type'] == 'tabular'
    assert actual['interpolation'] == interpolation
    np.testing.assert_array_equal(actual['x'], values[0])
    # Native Tabular normalizes by the supplied terminal CDF in its constructor.
    norm = values[2, -1]
    np.testing.assert_array_equal(actual['pdf'], values[1] / norm)
    np.testing.assert_array_equal(actual['cdf'], values[2] / norm)
    assert actual['integral'] == norm


def assert_law(actual, group):
    kind = decoded(group.attrs['type'])
    expected_type = {'continuous': 'continuous-tabular',
                     'discrete_photon': 'discrete-photon',
                     'level': 'level-inelastic'}.get(kind, kind)
    assert actual['type'] == expected_type
    if kind == 'uncorrelated':
        if 'energy' in group:
            assert_law(actual['energy'], group['energy'])
        else:
            assert actual['energy'] is None
        if 'angle' not in group:
            assert actual['angle'] == {'type': 'isotropic'}
        else:
            angle = group['angle']
            np.testing.assert_array_equal(actual['angle']['incident_energy'],
                                          angle['energy'][()])
            mu = angle['mu']
            offsets = list(mu.attrs['offsets']) + [mu.shape[1]]
            for i, table in enumerate(actual['angle']['tables']):
                assert_mu(table, mu[:, offsets[i]:offsets[i + 1]],
                          mu.attrs['interpolation'][i])
    elif kind in ('continuous', 'correlated', 'kalbach-mann'):
        np.testing.assert_array_equal(actual['incident_energy'], group['energy'][()])
        interp = group['energy'].attrs['interpolation']
        np.testing.assert_array_equal(actual['breakpoints'], interp[0])
        # int2interp strips the HDF5 two-dimensional interpolation qualifier.
        np.testing.assert_array_equal(actual['interpolation'], interp[1] % 10)
        data = group['energy_out' if kind == 'correlated' else 'distribution']
        offsets = list(data.attrs['offsets']) + [data.shape[1]]
        assert len(actual['tables']) == len(offsets) - 1
        for i, table in enumerate(actual['tables']):
            start, stop = offsets[i:i + 2]
            for j, key in enumerate(('energy', 'pdf', 'cdf')):
                np.testing.assert_array_equal(table[key], data[j, start:stop])
            assert table['n_discrete'] == data.attrs['n_discrete_lines'][i]
            assert table['interpolation'] == data.attrs['interpolation'][i]
            if kind == 'kalbach-mann':
                np.testing.assert_array_equal(table['r'], data[3, start:stop])
                np.testing.assert_array_equal(table['a'], data[4, start:stop])
            elif kind == 'correlated':
                assert len(table['angle']) == stop - start
                mu = group['mu']
                for j, angle in enumerate(table['angle'], start):
                    lo = int(data[4, j])
                    hi = int(data[4, j + 1]) if j + 1 < data.shape[1] else mu.shape[1]
                    assert_mu(angle, mu[:, lo:hi], max(1, int(data[3, j])))
    elif kind in ('maxwell', 'evaporation', 'watt'):
        assert actual['u'] == group.attrs['u']
        for name in ('a', 'b') if kind == 'watt' else ('theta',):
            assert_function(actual[name], group[name])
    else:
        names = {'discrete_photon': ('primary_flag', 'energy', 'atomic_weight_ratio'),
                 'level': ('threshold', 'mass_ratio'),
                 'nbody': ('n_particles', 'total_mass',
                           'atomic_weight_ratio', 'q_value')}
        for name in names[kind]:
            assert actual[name] == group.attrs[name]


@pytest.fixture(scope='module')
def native_data(tmp_path_factory, mpi_intracomm):
    library = openmc.data.DataLibrary.from_xml()
    names = [n for n in NUCLIDES if library.get_by_material(n) is not None]
    material = openmc.Material()
    for name in names:
        material.add_nuclide(name, 1.0)
    material.set_density('g/cm3', 1.0)
    sphere = openmc.Sphere(r=1, boundary_type='vacuum')
    cell = openmc.Cell(fill=material, region=-sphere)
    model = openmc.Model(geometry=openmc.Geometry([cell]))
    model.settings.run_mode = 'fixed source'
    model.settings.particles = 10
    model.settings.batches = 2
    model.settings.source = openmc.IndependentSource(space=openmc.stats.Point())
    directory = tmp_path_factory.mktemp('secondary')
    model.settings.output = {'path': str(directory), 'summary': False}
    try:
        model.init_lib(directory=str(directory),
                       output=False, intracomm=mpi_intracomm)
        yield library, names
    finally:
        model.finalize_lib()


@pytest.mark.parametrize('nuclide', NUCLIDES)
def test_discovery_metadata(native_data, nuclide):
    library, names = native_data
    if nuclide not in names:
        pytest.skip(f'{nuclide} not installed')
    products = openmc.lib.get_reaction_products(nuclide)
    assert products
    with h5py.File(library.get_by_material(nuclide)['path']) as f:
        reactions = f[nuclide]['reactions']
        expected_count = sum(sum(k.startswith('product_') for k in rx)
                             for rx in reactions.values())
        assert len(products) == expected_count
        for product in products:
            mt = product['mt']
            rx = reactions[f'reaction_{mt:03}']
            keys = sorted(k for k in rx if k.startswith('product_'))
            group = rx[keys[product['product_index']]]
            assert product['nuclide'] == nuclide
            assert product['schema_version'] == 1
            assert product['particle'] == decoded(group.attrs['particle'])
            assert product['emission_mode'] == decoded(group.attrs['emission_mode'])
            assert product['reference_frame'] == (
                'center-of-mass' if rx.attrs['center_of_mass'] else 'laboratory')
            assert_function(product['yield'], group['yield'])
            nd = group.attrs['n_distribution']
            assert len(product['distributions']) == nd
            assert len(product['applicability']) == (nd if nd > 1 else 0)
            for desc in product['distributions']:
                i = desc['distribution_index']
                d = group[f'distribution_{i}']
                assert desc['angle_energy_type'] == decoded(d.attrs['type'])
                assert desc['read_supported']
                assert not desc['mutation_supported']
                assert desc['mutation_reason']
                assert 'data' not in desc
                if nd > 1:
                    assert_function(product['applicability'][i], d['applicability'])
                assert desc['energy_mutation_supported'] == (
                    desc['energy_type'] == 'continuous-tabular')


@pytest.mark.parametrize('nuclide', NUCLIDES)
def test_complete_native_payloads(native_data, nuclide):
    library, names = native_data
    if nuclide not in names:
        pytest.skip(f'{nuclide} not installed')
    # One real example of every representation/particle/mode combination per
    # nuclide. Compare every nested outgoing table in each selected law.
    seen = set()
    with h5py.File(library.get_by_material(nuclide)['path']) as f:
        for p in openmc.lib.get_reaction_products(nuclide):
            for d in p['distributions']:
                key = (d['angle_energy_type'], d['energy_type'],
                       p['particle'], p['emission_mode'])
                if key in seen:
                    continue
                seen.add(key)
                actual = openmc.lib.get_emission_distribution(
                    nuclide, p['mt'], p['product_index'], d['distribution_index'])
                rx = f[nuclide]['reactions'][f"reaction_{p['mt']:03}"]
                keys = sorted(k for k in rx if k.startswith('product_'))
                product = rx[keys[p['product_index']]]
                group = product[f"distribution_{d['distribution_index']}"]
                assert_law(actual, group)


def test_copies_and_legacy_api(native_data):
    index = openmc.lib.nuclides['U235']._index
    new = openmc.lib.get_emission_distribution(index, 18)
    old = openmc.lib.get_emission_spectrum(index, 18)
    if new['energy']['type'] != 'continuous-tabular':
        pytest.skip('Installed U235 fission does not use ContinuousTabular')
    for before, after in zip(old['tables'], new['energy']['tables']):
        for key in ('energy', 'pdf', 'cdf'):
            np.testing.assert_array_equal(before[key], after[key])
    saved = deepcopy(new['energy']['tables'][0]['pdf'])
    new['energy']['tables'][0]['pdf'][:] = -123
    reread = openmc.lib.get_emission_distribution(index, 18)
    np.testing.assert_array_equal(reread['energy']['tables'][0]['pdf'], saved)
    product = openmc.lib.get_reaction_products('U235', 18, include_data=True)[0]
    np.testing.assert_array_equal(product['distributions'][0]['data']['energy']
                                  ['tables'][0]['pdf'], saved)


@pytest.mark.parametrize('args', [(-1, 18), (2**40, 18), ('U235', -1),
                                  ('U235', 18, -1), ('U235', 18, 0, -1)])
def test_invalid_python_selection(native_data, args):
    with pytest.raises(ValueError):
        openmc.lib.get_emission_distribution(*args)


@pytest.mark.parametrize('args', [(2**30, 18), ('U235', 999),
                                  ('U235', 18, 10000), ('U235', 18, 0, 10000)])
def test_invalid_native_selection(native_data, args):
    with pytest.raises(exceptions.OutOfBoundsError):
        openmc.lib.get_emission_distribution(*args)


def test_general_reader_observes_legacy_native_mutation(native_data):
    index = openmc.lib.nuclides['U235']._index
    original = openmc.lib.get_emission_spectrum(index, 18)
    changed = deepcopy(original)
    for table in changed['tables']:
        table['pdf'] *= np.exp(-table['energy'] / 3.e5)
    openmc.lib.simulation_init()
    try:
        openmc.lib.set_emission_spectrum(index, 18, changed)
        legacy = openmc.lib.get_emission_spectrum(index, 18)
        general = openmc.lib.get_emission_distribution(index, 18)['energy']
        assert not np.array_equal(original['tables'][0]['pdf'],
                                  general['tables'][0]['pdf'])
        for a, b in zip(legacy['tables'], general['tables']):
            for field in ('energy', 'pdf', 'cdf'):
                np.testing.assert_array_equal(a[field], b[field])
    finally:
        openmc.lib.set_emission_spectrum(index, 18, original)
        openmc.lib.simulation_finalize()
