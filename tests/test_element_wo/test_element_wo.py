#!/usr/bin/env python

import os
import sys

import numpy as np

sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, os.pardir))
from openmc import Material
from openmc.data import NATURAL_ABUNDANCE, atomic_mass


if __name__ == '__main__':
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
            assert np.isclose(densities[nuc][1], val, rtol=1.e-8)
        if nuc == 'O16':
            val = (NATURAL_ABUNDANCE[nuc] + NATURAL_ABUNDANCE['O18']) \
                  * atomic_mass(nuc) / water_am
            assert np.isclose(densities[nuc][1], val, rtol=1.e-8)
        if nuc == 'O17':
            val = NATURAL_ABUNDANCE[nuc] * atomic_mass(nuc) / water_am
            assert np.isclose(densities[nuc][1], val, rtol=1.e-8)
