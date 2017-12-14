#!/usr/bin/env python

from collections import Mapping
import os

import numpy as np
import pandas as pd
import pytest
import openmc.data


@pytest.fixture
def pu239():
    directory = os.path.dirname(os.environ['OPENMC_CROSS_SECTIONS'])
    filename = os.path.join(directory, 'Pu239.h5')
    return openmc.data.IncidentNeutron.from_hdf5(filename)


@pytest.fixture
def gd154():
    directory = os.environ['OPENMC_ENDF_DATA']
    filename = os.path.join(directory, 'neutrons', 'n-064_Gd_154.endf')
    return openmc.data.IncidentNeutron.from_endf(filename)


@pytest.fixture
def ace_file():
    directory = os.environ['OPENMC_ENDF_DATA']
    h1 = os.path.join(directory, 'neutrons', 'n-001_H_001.endf')
    retcode = openmc.data.njoy.make_ace(h1, ace='h1.ace')
    assert retcode == 0
    return os.path.join(os.getcwd(), 'h1.ace')


@pytest.fixture
def h1(ace_file):
    return openmc.data.IncidentNeutron.from_ace(ace_file)


def test_attributes(pu239):
    assert pu239.name == 'Pu239'
    assert pu239.mass_number == 239
    assert pu239.metastable == 0
    assert pu239.atomic_symbol == 'Pu'
    assert pu239.atomic_weight_ratio == pytest.approx(236.9986)


def test_energy_grid(pu239):
    assert isinstance(pu239.energy, Mapping)
    for temp, grid in pu239.energy.items():
        assert temp.endswith('K')
        assert np.all(np.diff(grid) >= 0.0)


def test_elastic(pu239):
    elastic = pu239.reactions[2]
    assert elastic.center_of_mass
    assert elastic.q_value == 0.0
    assert elastic.mt == 2
    assert '0K' in elastic.xs
    assert '294K' in elastic.xs
    assert len(elastic.products) == 1
    p = elastic.products[0]
    assert isinstance(p, openmc.data.Product)
    assert p.particle == 'neutron'
    assert p.emission_mode == 'prompt'
    assert len(p.distribution) == 1
    d = p.distribution[0]
    assert isinstance(d, openmc.data.UncorrelatedAngleEnergy)
    assert isinstance(d.angle, openmc.data.AngleDistribution)
    assert d.energy is None
    assert p.yield_(0.0) == 1.0


def test_fission(pu239):
    fission = pu239.reactions[18]
    assert not fission.center_of_mass
    assert fission.q_value == pytest.approx(198902000.0)
    assert fission.mt == 18
    assert '294K' in fission.xs
    assert len(fission.products) == 8
    prompt = fission.products[0]
    assert prompt.particle == 'neutron'
    assert prompt.yield_(1.0e-5) == pytest.approx(2.874262)
    delayed = [p for p in fission.products if p.emission_mode == 'delayed']
    assert len(delayed) == 6
    assert all(d.particle == 'neutron' for d in delayed)
    assert sum(d.decay_rate for d in delayed) == pytest.approx(4.037212)
    assert sum(d.yield_(1.0) for d in delayed) == pytest.approx(0.00645)
    photon = fission.products[-1]
    assert photon.particle == 'photon'


def test_get_reaction_components(h1):
    assert h1.get_reaction_components(1) == [2, 102]
    assert h1.get_reaction_components(101) == [102]
    assert h1.get_reaction_components(102) == [102]
    assert h1.get_reaction_components(51) == []


def test_resonances(gd154):
    res = gd154.resonances
    assert isinstance(res, openmc.data.Resonances)
    assert len(res.ranges) == 2
    resolved, unresolved = res.ranges
    assert isinstance(resolved, openmc.data.ReichMoore)
    assert isinstance(unresolved, openmc.data.Unresolved)
    assert resolved.energy_min == pytest.approx(1e-5)
    assert resolved.energy_max == pytest.approx(2760.)
    assert resolved.target_spin == 0.0
    assert resolved.channel_radius[0](1.0) == pytest.approx(0.74)
    assert isinstance(resolved.parameters, pd.DataFrame)
    assert (resolved.parameters['L'] == 0).all()
    assert (resolved.parameters['J'] <= 0.5).all()
    assert (resolved.parameters['fissionWidthA'] == 0.0).all()


def test_reconstruct(gd154):
    elastic = gd154.reactions[2].xs['0K']
    assert isinstance(elastic, openmc.data.ResonancesWithBackground)
    assert elastic(0.0253) == pytest.approx(5.7228949796394524)
