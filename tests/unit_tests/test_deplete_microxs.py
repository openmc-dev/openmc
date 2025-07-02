"""Basic unit tests for openmc.deplete.IndependentOperator instantiation

Modifies and resets environment variable OPENMC_CROSS_SECTIONS
to a custom file with new depletion_chain node
"""

from os import remove
from pathlib import Path

import pytest
import openmc
from openmc.deplete import MicroXS, get_microxs_from_multigroup, Chain
from openmc.mgxs import GROUP_STRUCTURES
import numpy as np

ONE_GROUP_XS = Path(__file__).parents[1] / "micro_xs_simple.csv"


def test_from_array():
    nuclides = [
        'U234',
        'U235',
        'U238',
        'U236',
        'O16',
        'O17',
        'I135',
        'Xe135',
        'Xe136',
        'Cs135',
        'Gd157',
        'Gd156']
    reactions = ['fission', '(n,gamma)']
    # These values are placeholders and are not at all
    # physically meaningful.
    data = np.array([[0.1, 0.],
                     [0.1, 0.],
                     [0.9, 0.],
                     [0.4, 0.],
                     [0., 0.],
                     [0., 0.],
                     [0., 0.1],
                     [0., 0.9],
                     [0., 0.],
                     [0., 0.],
                     [0., 0.1],
                     [0., 0.1]])
    data.shape = (12, 2, 1)

    MicroXS(data, nuclides, reactions)
    with pytest.raises(ValueError, match='Data array must be 3D'):
        MicroXS(data[:, 0], nuclides, reactions)


def test_csv():
    ref_xs = MicroXS.from_csv(ONE_GROUP_XS)
    ref_xs.to_csv('temp_xs.csv')
    temp_xs = MicroXS.from_csv('temp_xs.csv')
    assert np.all(ref_xs.data == temp_xs.data)
    remove('temp_xs.csv')


def test_from_multigroup_flux():
    energies = [0., 6.25e-1, 5.53e3, 8.21e5, 2.e7]
    flux = [1.1e-7, 1.2e-6, 1.3e-5, 1.4e-4]
    chain_file = Path(__file__).parents[1] / 'chain_simple.xml'
    kwargs = {'multigroup_flux': flux, 'chain_file': chain_file}

    # test with energy group structure from string
    microxs = MicroXS.from_multigroup_flux(energies='CASMO-4', **kwargs)
    assert isinstance(microxs, MicroXS)

    # test with energy group structure as floats
    microxs = MicroXS.from_multigroup_flux(energies=energies, **kwargs)
    assert isinstance(microxs, MicroXS)

    # test with nuclides provided
    microxs = MicroXS.from_multigroup_flux(
        energies=energies, nuclides=['Gd157', 'H1'], **kwargs
    )
    assert isinstance(microxs, MicroXS)
    assert microxs.nuclides == ['Gd157', 'H1']

    # test with reactions provided
    microxs = MicroXS.from_multigroup_flux(
        energies=energies, reactions=['fission', '(n,2n)'], **kwargs
    )
    assert isinstance(microxs, MicroXS)
    assert microxs.reactions == ['fission', '(n,2n)']


def test_multigroup_flux_same():
    chain_file = Path(__file__).parents[1] / 'chain_simple.xml'

    # Generate micro XS based on 4-group flux
    energies = [0., 6.25e-1, 5.53e3, 8.21e5, 2.e7]
    flux_per_ev = [0.3, 0.3, 1.0, 1.0]
    flux = flux_per_ev * np.diff(energies)
    flux_sum = flux.sum()
    microxs_4g = MicroXS.from_multigroup_flux(
        energies=energies, multigroup_flux=flux, chain_file=chain_file)

    # from_multigroup_flux should not modify the flux
    assert flux.sum() == flux_sum

    # Generate micro XS based on 2-group flux, where the boundaries line up with
    # the 4 group flux and have the same flux per eV across the full energy
    # range
    energies = [0., 5.53e3, 2.0e7]
    flux_per_ev = [0.3, 1.0]
    flux = flux_per_ev * np.diff(energies)
    microxs_2g = MicroXS.from_multigroup_flux(
        energies=energies, multigroup_flux=flux, chain_file=chain_file)

    assert microxs_4g.data == pytest.approx(microxs_2g.data)