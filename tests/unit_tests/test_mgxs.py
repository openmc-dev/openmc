import numpy as np
import pytest

import openmc


@pytest.mark.parametrize('mgxs_class', [
    openmc.mgxs.PhotonTransferMatrixXS,
    openmc.mgxs.NuFissionMatrixXS,
    openmc.mgxs.ScatterMatrixXS,
])
@pytest.mark.parametrize('in_groups, out_groups, expected', [
    ([1], [2, 3], [[8.0, 7.0]]),
    ([1], [], [[9.0, 8.0, 7.0]]),
    ([], [2, 3], [[8.0, 7.0], [5.0, 4.0], [2.0, 1.0]]),
    ([1, 3], [1, 3], [[9.0, 7.0], [3.0, 1.0]]),
])
def test_matrix_slice_group_selection(mgxs_class, in_groups, out_groups,
                                     expected):
    groups = openmc.mgxs.EnergyGroups([1.0, 10.0, 100.0, 1000.0])
    mgxs = mgxs_class(openmc.Material(), 'material', groups)
    if isinstance(mgxs, openmc.mgxs.ScatterMatrixXS):
        mgxs.correction = None
    for key, tally in mgxs.tallies.items():
        tally._derived = True
        tally._mean = np.zeros(tally.shape)
        tally._std_dev = np.zeros(tally.shape)
        if key == 'flux':
            tally._mean.fill(1.0)
    key = ('primary photon production'
           if isinstance(mgxs, openmc.mgxs.PhotonTransferMatrixXS)
           else mgxs.rxn_type)
    mgxs.tallies[key]._mean[:, 0, 0] = np.arange(1.0, 10.0)

    sliced = mgxs.get_slice(in_groups=in_groups, out_groups=out_groups)
    expected = np.array(expected)
    assert np.array_equal(sliced.get_xs(), expected)
    assert np.array_equal(sliced.get_xs(row_column='outin'), expected.T)
    assert np.array_equal(
        sliced.get_xs(order_groups='decreasing'), expected[::-1, ::-1])
    assert np.array_equal(
        sliced.get_xs(in_groups=[1], out_groups=[expected.shape[1]]),
        expected[:1, -1:])

    resliced = sliced.get_slice(
        in_groups=[expected.shape[0]], out_groups=[1])
    assert np.array_equal(resliced.get_xs(), expected[-1:, :1])
    assert np.array_equal(
        resliced.get_xs(in_groups=[1], out_groups=[1]), expected[-1:, :1])
    assert np.array_equal(
        mgxs.get_xs(), [[9.0, 8.0, 7.0], [6.0, 5.0, 4.0],
                        [3.0, 2.0, 1.0]])

    for kwargs in ({'in_groups': [expected.shape[0] + 1]},
                   {'out_groups': [expected.shape[1] + 1]},
                   {'in_groups': [0]}, {'out_groups': [0]}):
        with pytest.raises(ValueError):
            sliced.get_xs(**kwargs)
        with pytest.raises(ValueError):
            sliced.get_slice(**kwargs)


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
