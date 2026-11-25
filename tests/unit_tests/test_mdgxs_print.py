"""Reports for delayed cross sections retain the correct data dimensions."""

import re

import numpy as np
import pytest

import openmc
import openmc.mgxs


@pytest.mark.parametrize('mgxs_type', openmc.mgxs.MDGXS_TYPES)
@pytest.mark.parametrize('delayed_groups', [None, [3], [1, 3]])
@pytest.mark.parametrize('angles', [(1, 1), (1, 2), (2, 1), (2, 2)])
@pytest.mark.parametrize('num_groups', [1, 2])
@pytest.mark.parametrize('nuclides', [None, 'all', 'sum', ['U238']])
def test_mdgxs_print_xs(mgxs_type, delayed_groups, angles, num_groups,
                        nuclides, capsys, xs_type='macro'):
    material = openmc.Material()
    material.add_nuclide('U235', 0.25)
    material.add_nuclide('U238', 0.75)
    material.set_density('atom/b-cm', 0.02)
    energy_groups = openmc.mgxs.EnergyGroups(
        np.linspace(0.0, 20.0e6, num_groups + 1))
    mgxs = openmc.mgxs.MDGXS.get_mgxs(
        mgxs_type, material, 'material', energy_groups, delayed_groups,
        by_nuclide=nuclides is not None,
        num_polar=angles[0], num_azimuthal=angles[1])

    # Exercise real tally arithmetic and get_xs, without a transport run.
    for i, tally in enumerate(mgxs.tallies.values(), 1):
        tally._derived = True
        tally._mean = np.arange(1, np.prod(tally.shape) + 1,
                                dtype=float).reshape(tally.shape) + i
        tally._std_dev = tally._mean * (i * 0.01)

    expected = []
    if nuclides == 'all':
        query_nuclides = ['U235', 'U238']
    elif nuclides is None or nuclides == 'sum':
        query_nuclides = ['sum']
    else:
        query_nuclides = nuclides
    for nuclide in query_nuclides:
        selection = 'sum' if nuclide == 'sum' else [nuclide]
        for group in delayed_groups or [None]:
            selected_groups = 'all' if group is None else [group]
            means = mgxs.get_xs(nuclides=selection,
                                delayed_groups=selected_groups, value='mean',
                                xs_type=xs_type)
            errors = mgxs.get_xs(nuclides=selection,
                                 delayed_groups=selected_groups,
                                 value='rel_err', xs_type=xs_type)
            assert np.all(np.isfinite(means))
            assert np.all(means < 1.0e6)
            expected.extend(f'{mean:.2e} +/- {100.0 * error:.2e}%'
                            for mean, error in zip(means.flat, errors.flat))

    mgxs.print_xs(nuclides='all' if nuclides is None else nuclides,
                  xs_type=xs_type)
    output = capsys.readouterr().out
    values = re.findall(r'[\d.e+-]+ \+/- [\d.e+-]+%', output)
    assert values == expected
    for group in delayed_groups or []:
        assert f'Delayed Group {group}:' in output
    if delayed_groups is None:
        assert 'Delayed Group' not in output
    if mgxs_type == 'decay-rate':
        assert 'eV]' not in output
        assert 'Decay Rates [s^-1]' in output
    else:
        assert 'Group 1' in output
        assert f'[{mgxs.get_units(xs_type)}]' in output


@pytest.mark.parametrize('mgxs_type', ['delayed-nu-fission',
                                       'delayed-nu-fission matrix'])
@pytest.mark.parametrize('nuclides', ['all', 'sum'])
def test_mdgxs_print_micro(mgxs_type, nuclides, capsys):
    test_mdgxs_print_xs(mgxs_type, [1, 3], (2, 2), 2, nuclides, capsys,
                        xs_type='micro')


@pytest.mark.parametrize('mgxs_type', openmc.mgxs.MDGXS_TYPES)
def test_mdgxs_print_mesh(mgxs_type, capsys):
    mesh = openmc.RegularMesh()
    mesh.dimension = (2,)
    mesh.lower_left = (0.0,)
    mesh.upper_right = (2.0,)
    mgxs = openmc.mgxs.MDGXS.get_mgxs(
        mgxs_type, mesh, 'mesh', openmc.mgxs.EnergyGroups([0.0, 1.0, 20.0e6]),
        delayed_groups=[1, 3])
    for i, tally in enumerate(mgxs.tallies.values(), 1):
        tally._derived = True
        tally._mean = np.arange(1, np.prod(tally.shape) + 1,
                                dtype=float).reshape(tally.shape) + i
        tally._std_dev = tally._mean * (i * 0.01)

    expected = []
    for subdomain in [(1,), (2,)]:
        for group in [1, 3]:
            means = mgxs.get_xs(subdomains=[subdomain], delayed_groups=[group])
            errors = mgxs.get_xs(subdomains=[subdomain],
                                 delayed_groups=[group],
                                 value='rel_err')
            expected.extend(f'{mean:.2e} +/- {100.0 * error:.2e}%'
                            for mean, error in zip(means.flat, errors.flat))
    mgxs.print_xs()
    output = capsys.readouterr().out
    assert re.findall(r'[\d.e+-]+ \+/- [\d.e+-]+%', output) == expected
    assert output.count('Delayed Group 1:') == 2
    assert output.count('Delayed Group 3:') == 2


@pytest.mark.parametrize('delayed_groups', [None, [1, 3]])
@pytest.mark.parametrize('angles', [(1, 1), (2, 2)])
def test_decay_rate_print_xs_type(delayed_groups, angles, capsys):
    material = openmc.Material()
    material.add_nuclide('U235', 0.25)
    material.add_nuclide('U238', 0.75)
    material.set_density('atom/b-cm', 0.02)
    mgxs = openmc.mgxs.DecayRate(
        material, 'material', openmc.mgxs.EnergyGroups([0.0, 1.0, 20.0e6]),
        delayed_groups=delayed_groups, by_nuclide=True,
        num_polar=angles[0], num_azimuthal=angles[1])
    for i, tally in enumerate(mgxs.tallies.values(), 1):
        tally._derived = True
        tally._mean = np.arange(1, np.prod(tally.shape) + 1,
                                dtype=float).reshape(tally.shape) + i
        tally._std_dev = tally._mean * (i * 0.01)

    mgxs.print_xs(nuclides='sum', xs_type='macro')
    macro_output = capsys.readouterr().out
    mgxs.print_xs(nuclides='sum', xs_type='micro')
    micro_output = capsys.readouterr().out
    assert macro_output == micro_output
    assert 'Decay Rates [s^-1]' in macro_output
    assert 'eV]' not in macro_output
    expected_bins = len(delayed_groups or [None]) * angles[0] * angles[1]
    assert macro_output.count('Decay Rate:') == expected_bins
