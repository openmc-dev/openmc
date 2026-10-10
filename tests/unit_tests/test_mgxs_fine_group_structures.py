"""Tests for openmc.mgxs.GROUP_STRUCTURES and build_fine_group_structure."""

from pathlib import Path

import numpy as np
import pytest

import openmc
from openmc.mgxs import GROUP_STRUCTURES, build_fine_group_structure

# Edges tabulated from: !!TODOXYZ!!
REFERENCE_FILE = Path(__file__).with_name('fine_group_structures_reference.npz')

# VESTA-43000 macro group bounds [eV] and number of fine group bins per macro group
VESTA43000_MACRO_BOUNDS = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1., 10., 100., 1e3, 1e4, 1e5,
                1e6, 1e7, 2e7]
VESTA43000_NUM_FINE_BINS = [1000, 1000, 1000, 1000, 1000, 4000, 4000, 10000,
                    10000, 4000, 4000, 1000, 1000]


@pytest.fixture(scope='module')
def reference():
    with np.load(REFERENCE_FILE) as data:
        return {name: data[name] for name in data.files}


@pytest.mark.parametrize('name', sorted(GROUP_STRUCTURES))
def test_group_structure_well_formed(name):
    edges = GROUP_STRUCTURES[name]
    assert isinstance(edges, np.ndarray)
    assert edges.ndim == 1
    assert edges.dtype == np.float64
    assert np.all(np.isfinite(edges))
    assert np.all(np.diff(edges) > 0.0)
    assert edges[0] >= 0.0


@pytest.mark.parametrize('spacing', ['log', 'linear'])
def test_build_scalar_arguments(spacing):
    macro_bounds = [1e-5, 1., 1e3, 2e7]
    expected = build_fine_group_structure(macro_bounds, [10, 10, 10], [spacing] * 3)
    edges = build_fine_group_structure(macro_bounds, 10, spacing)
    np.testing.assert_array_equal(edges, expected)


def test_build_mixed_spacing():
    macro_bounds = [0., 1., 1e3, 2e7]
    num_fine_bins = [4, 30, 7]
    spacing = ['linear', 'log', 'linear']
    edges = build_fine_group_structure(macro_bounds, num_fine_bins, spacing)
    assert edges.size == sum(num_fine_bins) + 1

    # Macro group bounds are reproduced exactly
    idx = np.cumsum([0] + num_fine_bins)
    np.testing.assert_array_equal(edges[idx], macro_bounds)

    # Constant lethargy ('log') or energy ('linear') width within a macro group
    for kind, start, stop in zip(spacing, idx[:-1], idx[1:]):
        segment = edges[start:stop + 1]
        if kind == 'log':
            widths = np.diff(np.log(segment))
        else:
            widths = np.diff(segment)
        np.testing.assert_allclose(widths, widths[0], rtol=1e-9)


def test_build_linear_from_zero():
    edges = build_fine_group_structure([0., 1.], 4, 'linear')
    np.testing.assert_array_equal(edges, [0., 0.25, 0.5, 0.75, 1.])


def test_build_uint8_num_fine_bins():
    # np.uint8(255) + 1 would wrap to 0 if the count stayed a numpy scalar
    edges = build_fine_group_structure([1., 10.], np.uint8(255))
    assert edges.size == 256
    assert edges[0] == 1.
    assert edges[-1] == 10.


def test_build_array_num_fine_bins_tuple_spacing():
    edges = build_fine_group_structure(
        [1., 10., 100.], np.array([3, 4], dtype=np.int64),
        spacing=('log', 'linear'))
    assert edges.size == 8
    np.testing.assert_array_equal(edges[[0, 3, 7]], [1., 10., 100.])


@pytest.mark.parametrize('macro_bounds, num_fine_bins, spacing', [
    pytest.param([1., 10., 100.], [10], 'log', id='num-groups-length'),
    pytest.param([1., 10., 100.], 10, ['log'], id='spacing-length'),
    pytest.param([1., 100., 10.], 10, 'log', id='decreasing-bounds'),
    pytest.param([1., 1., 10.], 10, 'log', id='repeated-bound'),
    pytest.param([1.], 10, 'log', id='single-bound'),
    pytest.param([[1., 10.]], 10, 'log', id='2d-bounds'),
    pytest.param([-1., 1.], 2, 'linear', id='negative-bound'),
    pytest.param([1., np.inf], 2, 'log', id='infinite-bound'),
    pytest.param([1., 10.], 0, 'log', id='zero-groups'),
    pytest.param([1., 10.], 10, 'cubic', id='unknown-spacing'),
    pytest.param([0., 10.], 10, 'log', id='log-from-zero'),
    pytest.param([1., 1. + 1e-12], 100000, 'log', id='too-narrow-macrobound'),
])
def test_build_invalid(macro_bounds, num_fine_bins, spacing):
    with pytest.raises(ValueError):
        build_fine_group_structure(macro_bounds, num_fine_bins, spacing)


def test_build_non_integer_num_fine_bins():
    with pytest.raises(TypeError):
        build_fine_group_structure([1., 10.], 2.5)


@pytest.mark.parametrize('name, rtol', [
    ('VESTA-43000', 1e-12),
    # FOMG edges were tabulated to 7 significant figures
    ('FOMG-16000', 1e-6),
])
def test_reference_edges(reference, name, rtol):
    edges = GROUP_STRUCTURES[name]
    expected = reference[name]
    assert edges.size == expected.size
    assert edges[0] == expected[0]
    assert edges[-1] == expected[-1]
    np.testing.assert_allclose(edges, expected, rtol=rtol, atol=0)


def test_vesta_43000():
    edges = GROUP_STRUCTURES['VESTA-43000']
    assert edges.size == 43001

    # Macrogroup endpoints are exact; lethargy width ln(hi/lo)/n inside
    idx = np.cumsum([0] + VESTA43000_NUM_FINE_BINS)
    np.testing.assert_array_equal(edges[idx], VESTA43000_MACRO_BOUNDS)
    for lo, hi, n, start, stop in zip(VESTA43000_MACRO_BOUNDS[:-1], VESTA43000_MACRO_BOUNDS[1:],
                                      VESTA43000_NUM_FINE_BINS, idx[:-1], idx[1:]):
        du = np.diff(np.log(edges[start:stop + 1]))
        np.testing.assert_allclose(du, np.log(hi / lo) / n, rtol=1e-9)


def test_fomg_16000():
    edges = GROUP_STRUCTURES['FOMG-16000']
    assert edges.size == 16001
    assert edges[0] == 1e-5
    assert edges[1000] == 1.0
    assert edges[15000] == 2e6
    assert edges[-1] == 1.96e7

    # 1 meV bins from the second edge up to 1 eV
    np.testing.assert_allclose(np.diff(edges[1:1001]), 1e-3, rtol=1e-9)
    # Equal-lethargy bins from 1 eV to 2 MeV
    np.testing.assert_allclose(np.diff(np.log(edges[1000:15001])),
                               np.log(2e6) / 14000, rtol=1e-9)
    # 17.6 keV bins from 2 MeV to 19.6 MeV
    np.testing.assert_allclose(np.diff(edges[15000:]), 17.6e3, rtol=1e-9)


def test_vesta_100000():
    # Check the definition only
    edges = GROUP_STRUCTURES['VESTA-100000']
    assert edges.size == 100001
    assert edges[0] == 1e-5
    assert edges[-1] == 2e7
    np.testing.assert_allclose(np.diff(np.log(edges)),
                               np.log(2e7 / 1e-5) / 100000, rtol=1e-9)


@pytest.mark.parametrize('name, n_groups', [
    ('FOMG-16000', 16000),
    ('VESTA-43000', 43000),
    ('VESTA-100000', 100000),
])
def test_name_lookup(name, n_groups):
    assert openmc.EnergyFilter.from_group_structure(name).num_bins == n_groups
    assert openmc.mgxs.EnergyGroups(name).num_fine_bins == n_groups
    f = openmc.ParticleProductionFilter('photon', energies=name)
    assert f.num_energy_bins == n_groups
