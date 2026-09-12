import numpy as np

import openmc


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
