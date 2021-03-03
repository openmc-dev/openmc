#!/usr/bin/env python

import os
import sys

import pytest

from openmc import Material
from openmc.data import NATURAL_ABUNDANCE, atomic_mass


def test_element_wo():
    # This test doesn't require an OpenMC run.  We just need to make sure the
    # element.expand() method expands elements with the proper nuclide
    # compositions.

    h_am = (NATURAL_ABUNDANCE['H1'] * atomic_mass('H1') +
            NATURAL_ABUNDANCE['H2'] * atomic_mass('H2'))
    o_am = (NATURAL_ABUNDANCE['O17'] * atomic_mass('O17') +
            (NATURAL_ABUNDANCE['O16'] + NATURAL_ABUNDANCE['O18'])
            * atomic_mass('O16'))
    water_am = 2 * h_am + o_am

    water = Material()
    water.add_element('O', o_am / water_am, 'wo')
    water.add_element('H', 2 * h_am / water_am, 'wo')
    densities = water.get_nuclide_densities()

    for nuc in densities.keys():
        assert nuc in ('H1', 'H2', 'O16', 'O17')

        if nuc in ('H1', 'H2'):
            val = 2 * NATURAL_ABUNDANCE[nuc] * atomic_mass(nuc) / water_am
            assert densities[nuc][1] == pytest.approx(val)
        if nuc == 'O16':
            val = (NATURAL_ABUNDANCE[nuc] + NATURAL_ABUNDANCE['O18']) \
                  * atomic_mass(nuc) / water_am
            assert densities[nuc][1] == pytest.approx(val)
        if nuc == 'O17':
            val = NATURAL_ABUNDANCE[nuc] * atomic_mass(nuc) / water_am
            assert densities[nuc][1] == pytest.approx(val)
