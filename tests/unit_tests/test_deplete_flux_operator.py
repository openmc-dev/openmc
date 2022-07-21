"""Basic unit tests for openmc.deplete.FluxDepletionOperator instantiation

Modifies and resets environment variable OPENMC_CROSS_SECTIONS
to a custom file with new depletion_chain node
"""

from pathlib import Path

import pytest
from openmc.deplete.flux_operator import FluxDepletionOperator
import pandas as pd
import numpy as np

FLUX = 5e16
CHAIN_PATH = Path(__file__).parents[1] / "chain_simple.xml"
ONE_GROUP_XS = Path(__file__).parents[1] / "micro_xs_simple.csv"


def test_create_micro_xs_from_data_array():
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

    FluxDepletionOperator.create_micro_xs_from_data_array(
        nuclides, reactions, data)
    with pytest.raises(ValueError, match=r'Nuclides list of length \d* and '
                       r'reactions array of length \d* do not '
                       r'match dimensions of data array of shape \(\d*\,d*\)'):
        FluxDepletionOperator.create_micro_xs_from_data_array(
            nuclides, reactions, data[:, 0])


def test_create_micro_xs_from_csv():
    FluxDepletionOperator.create_micro_xs_from_csv(ONE_GROUP_XS)


def test_operator_init():
    """The test uses a temporary dummy chain. This file will be removed
    at the end of the test, and only contains a depletion_chain node."""
    nuclides = {'U234': 8.922411359424315e+18,
                'U235': 9.98240191860822e+20,
                'U238': 2.2192386373095893e+22,
                'U236': 4.5724195495061115e+18,
                'O16': 4.639065406771322e+22,
                'O17': 1.7588724018066158e+19}
    micro_xs = FluxDepletionOperator.create_micro_xs_from_csv(ONE_GROUP_XS)
    nuclide_flux_operator = FluxDepletionOperator.from_nuclides(
        1, nuclides, micro_xs, FLUX_SPECTRA, CHAIN_PATH)
