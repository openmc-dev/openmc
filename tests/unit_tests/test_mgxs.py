import numpy as np
import pytest

import openmc


@pytest.fixture(params=[
    'scatter matrix', 'consistent scatter matrix',
    'scatter probability matrix',
])
def scattering_matrix(request):
    groups = openmc.mgxs.EnergyGroups([1.0, 10.0, 100.0, 1000.0])
    mgxs = openmc.mgxs.MGXS.get_mgxs(
        request.param, openmc.Material(), 'material', groups)
    if request.param.startswith('consistent'):
        mgxs.correction = None
    for tally in mgxs.tallies.values():
        tally._derived = True
        tally._mean = np.ones(tally.shape)
        tally._std_dev = np.zeros(tally.shape)
    return mgxs


def test_matrix_slice_requires_full_groups(scattering_matrix):
    mgxs = scattering_matrix
    with pytest.raises(ValueError, match=r'Use get_xs\(\)'):
        mgxs.get_slice(out_groups=[1])
    assert np.allclose(
        mgxs.get_slice(in_groups=[3, 2, 1], out_groups=[3, 2, 1]).get_xs(),
        mgxs.get_xs())


@pytest.mark.parametrize('scattering_matrix', [
    'consistent scatter matrix', 'scatter probability matrix',
], indirect=True)
def test_matrix_slice_retains_normalization(scattering_matrix):
    full = scattering_matrix.get_xs()
    sliced = scattering_matrix.get_slice(in_groups=[1, 2])
    assert np.allclose(sliced.get_xs(), full[:2, :])
    assert np.allclose(
        sliced.get_slice(in_groups=[1], out_groups=[1, 2, 3]).get_xs(),
        full[:1, :])


@pytest.mark.parametrize('scattering_matrix', ['scatter matrix'],
                         indirect=True)
def test_matrix_slice_to_p0(scattering_matrix):
    mgxs = scattering_matrix
    mgxs.correction = None
    mgxs.legendre_order = 1
    with pytest.warns(UserWarning, match='ignored'):
        mgxs.correction = 'P0'
    with pytest.raises(ValueError, match=r'Use get_xs\(\)'):
        mgxs.get_slice(in_groups=[1], legendre_order=0)


def test_condensation_with_shared_energy_filter():
    groups = openmc.mgxs.EnergyGroups([1.0, 10.0, 100.0, 1000.0])
    total = openmc.mgxs.TotalXS(openmc.Material(), 'material', groups)
    flux = total.tallies['flux']
    reaction = total.tallies['total']
    energy = flux.find_filter(openmc.EnergyFilter)
    assert reaction.find_filter(openmc.EnergyFilter) is energy

    for tally, mean, std_dev in (
            (flux, [1.0, 2.0, 3.0], [0.1, 0.2, 0.3]),
            (reaction, [2.0, 8.0, 18.0], [0.2, 0.4, 0.6])):
        tally._derived = True
        tally._mean = np.array(mean).reshape(tally.shape)
        tally._std_dev = np.array(std_dev).reshape(tally.shape)

    coarse_groups = openmc.mgxs.EnergyGroups([1.0, 100.0, 1000.0])
    condensed = total.get_condensed_xs(coarse_groups)
    assert np.allclose(condensed.get_xs(), [6.0, 10.0 / 3.0])
    assert np.array_equal(
        condensed.tallies['flux'].mean.ravel(), [3.0, 3.0])
    assert np.allclose(
        condensed.tallies['total'].std_dev.ravel(),
        [np.hypot(0.2, 0.4), 0.6])
    for tally in condensed.tallies.values():
        assert np.array_equal(
            tally.find_filter(openmc.EnergyFilter).values,
            coarse_groups.group_edges)

    assert reaction.find_filter(openmc.EnergyFilter) is energy
    assert np.array_equal(energy.values, groups.group_edges)
    assert np.array_equal(total.get_xs(), [6.0, 4.0, 2.0])
    assert np.array_equal(flux.mean.ravel(), [1.0, 2.0, 3.0])

    # A second condensation must still use the original fine-group data.
    single_group = openmc.mgxs.EnergyGroups([1.0, 1000.0])
    assert np.allclose(total.get_condensed_xs(single_group).get_xs(), 28.0 / 6.0)
