#!/usr/bin/env python

from collections.abc import Mapping
import os
from math import log

import numpy as np
import pytest
from uncertainties import ufloat
import openmc.data


def ufloat_close(a, b):
    assert a.nominal_value == pytest.approx(b.nominal_value)
    assert a.std_dev == pytest.approx(b.std_dev)


@pytest.fixture(scope='module')
def nb90():
    """Nb90 decay data."""
    endf_data = os.environ['OPENMC_ENDF_DATA']
    filename = os.path.join(endf_data, 'decay', 'dec-041_Nb_090.endf')
    return openmc.data.Decay.from_endf(filename)


@pytest.fixture(scope='module')
def u235_yields():
    """U235 fission product yield data."""
    endf_data = os.environ['OPENMC_ENDF_DATA']
    filename = os.path.join(endf_data, 'nfy', 'nfy-092_U_235.endf')
    return openmc.data.FissionProductYields.from_endf(filename)


def test_nb90_halflife(nb90):
    ufloat_close(nb90.half_life, ufloat(52560.0, 180.0))
    ufloat_close(nb90.decay_constant, log(2.)/nb90.half_life)


def test_nb90_nuclide(nb90):
    assert nb90.nuclide['atomic_number'] == 41
    assert nb90.nuclide['mass_number'] == 90
    assert nb90.nuclide['isomeric_state'] == 0
    assert nb90.nuclide['parity'] == 1.0
    assert nb90.nuclide['spin'] == 8.0
    assert not nb90.nuclide['stable']
    assert nb90.nuclide['mass'] == pytest.approx(89.13888)


def test_nb90_modes(nb90):
    assert len(nb90.modes) == 2
    ec = nb90.modes[0]
    assert ec.modes == ['ec/beta+']
    assert ec.parent == 'Nb90'
    assert ec.daughter == 'Zr90'
    assert 'Nb90 -> Zr90' in str(ec)
    ufloat_close(ec.branching_ratio, ufloat(0.0147633, 0.0003386195))
    ufloat_close(ec.energy, ufloat(6111000., 4000.))

    # Make sure branching ratios sum to 1
    total = sum(m.branching_ratio for m in nb90.modes)
    assert total.nominal_value == pytest.approx(1.0)


def test_nb90_spectra(nb90):
    assert sorted(nb90.spectra.keys()) == ['e-', 'ec/beta+', 'gamma', 'xray']


def test_fpy(u235_yields):
    assert u235_yields.nuclide['atomic_number'] == 92
    assert u235_yields.nuclide['mass_number'] == 235
    assert u235_yields.nuclide['isomeric_state'] == 0
    assert u235_yields.nuclide['name'] == 'U235'
    assert u235_yields.energies == pytest.approx([0.0253, 500.e3, 1.4e7])

    assert len(u235_yields.cumulative) == 3
    thermal = u235_yields.cumulative[0]
    ufloat_close(thermal['I135'], ufloat(0.0628187, 0.000879461))

    assert len(u235_yields.independent) == 3
    thermal = u235_yields.independent[0]
    ufloat_close(thermal['I135'], ufloat(0.0292737, 0.000819663))
