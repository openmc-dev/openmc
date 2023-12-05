"""Basic unit tests for openmc.deplete.IndependentOperator instantiation

Modifies and resets environment variable OPENMC_CROSS_SECTIONS
to a custom file with new depletion_chain node
"""

from os import remove
from pathlib import Path

import pytest
from openmc.deplete import MicroXS
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
    with pytest.raises(ValueError, match=r'Nuclides list of length \d* and '
                       r'reactions array of length \d* do not '
                       r'match dimensions of data array of shape \(\d*\, \d*\)'):
        MicroXS(data[:, 0], nuclides, reactions)

def test_csv():
    ref_xs = MicroXS.from_csv(ONE_GROUP_XS)
    ref_xs.to_csv('temp_xs.csv')
    temp_xs = MicroXS.from_csv('temp_xs.csv')
    assert np.all(ref_xs.data == temp_xs.data)
    remove('temp_xs.csv')

def test_from_multi_group_flux():

    # test with energy group structure from string
    microxs = MicroXS.from_multi_group_flux(
        energies='CASMO-4',
        multi_group_flux=[1.1e-7, 1.2e-6, 1.3e-5, 1.4e-4],
        temperature=293,
        chain_file=Path(__file__).parents[1] / 'chain_simple.xml'
    )
    assert isinstance(microxs, MicroXS)

    # test with energy group structure as floats
    microxs = MicroXS.from_multi_group_flux(
        energies=[0., 6.25e-1, 5.53e3, 8.21e5, 2.e7],
        multi_group_flux=[1.1e-7, 1.2e-6, 1.3e-5, 1.4e-4],
        temperature=293,
        chain_file=Path(__file__).parents[1] / 'chain_simple.xml'
    )
    assert isinstance(microxs, MicroXS)

    # test with nuclides provided
    microxs = MicroXS.from_multi_group_flux(
        energies=[0., 6.25e-1, 5.53e3, 8.21e5, 2.e7],
        multi_group_flux=[1.1e-7, 1.2e-6, 1.3e-5, 1.4e-4],
        temperature=293,
        chain_file=Path(__file__).parents[1] / 'chain_simple.xml',
        nuclides=['Gd157']
    )
    assert microxs.nuclides == ['Gd157']
    assert isinstance(microxs, MicroXS)
