#!/usr/bin/env python

from collections import Mapping, Callable
import os

import numpy as np
import pandas as pd
import pytest
import openmc.data


_TEMPERATURES = [300., 600., 900.]
_ENDF_DATA = os.environ['OPENMC_ENDF_DATA']


@pytest.fixture(scope='module')
def pu239():
    """Pu239 HDF5 data."""
    directory = os.path.dirname(os.environ['OPENMC_CROSS_SECTIONS'])
    filename = os.path.join(directory, 'Pu239.h5')
    return openmc.data.IncidentNeutron.from_hdf5(filename)


@pytest.fixture(scope='module')
def xe135():
    """Xe135 ENDF data (contains SLBW resonance range)"""
    filename = os.path.join(_ENDF_DATA, 'neutrons', 'n-054_Xe_135.endf')
    return openmc.data.IncidentNeutron.from_endf(filename)


@pytest.fixture(scope='module')
def sm150():
    """Sm150 ENDF data (contains MLBW resonance range)"""
    filename = os.path.join(_ENDF_DATA, 'neutrons', 'n-062_Sm_150.endf')
    return openmc.data.IncidentNeutron.from_endf(filename)


@pytest.fixture(scope='module')
def gd154():
    """Gd154 ENDF data (contains Reich Moore resonance range)"""
    filename = os.path.join(_ENDF_DATA, 'neutrons', 'n-064_Gd_154.endf')
    return openmc.data.IncidentNeutron.from_endf(filename)


@pytest.fixture(scope='module')
def cl35():
    """Cl35 ENDF data (contains RML resonance range)"""
    filename = os.path.join(_ENDF_DATA, 'neutrons', 'n-017_Cl_035.endf')
    return openmc.data.IncidentNeutron.from_endf(filename)


@pytest.fixture(scope='module')
def am241():
    """Am241 ENDF data (contains Madland-Nix fission energy distribution)."""
    filename = os.path.join(_ENDF_DATA, 'neutrons', 'n-095_Am_241.endf')
    return openmc.data.IncidentNeutron.from_endf(filename)


@pytest.fixture(scope='module')
def u233():
    """U233 ENDF data (contains Watt fission energy distribution)."""
    filename = os.path.join(_ENDF_DATA, 'neutrons', 'n-092_U_233.endf')
    return openmc.data.IncidentNeutron.from_endf(filename)


@pytest.fixture(scope='module')
def u236():
    """U236 ENDF data (contains Watt fission energy distribution)."""
    filename = os.path.join(_ENDF_DATA, 'neutrons', 'n-092_U_236.endf')
    return openmc.data.IncidentNeutron.from_endf(filename)


@pytest.fixture(scope='module')
def na22():
    """Na22 ENDF data (contains evaporation spectrum)."""
    filename = os.path.join(_ENDF_DATA, 'neutrons', 'n-011_Na_022.endf')
    return openmc.data.IncidentNeutron.from_endf(filename)


@pytest.fixture(scope='module')
def be9():
    """Be9 ENDF data (contains laboratory angle-energy distribution)."""
    filename = os.path.join(_ENDF_DATA, 'neutrons', 'n-004_Be_009.endf')
    return openmc.data.IncidentNeutron.from_endf(filename)


@pytest.fixture(scope='module')
def h2():
    endf_file = os.path.join(_ENDF_DATA, 'neutrons', 'n-001_H_002.endf')
    return openmc.data.IncidentNeutron.from_njoy(
        endf_file, temperatures=_TEMPERATURES)


@pytest.fixture(scope='module')
def am244():
    endf_file = os.path.join(_ENDF_DATA, 'neutrons', 'n-095_Am_244.endf')
    return openmc.data.IncidentNeutron.from_njoy(endf_file)


def test_attributes(pu239):
    assert pu239.name == 'Pu239'
    assert pu239.mass_number == 239
    assert pu239.metastable == 0
    assert pu239.atomic_symbol == 'Pu'
    assert pu239.atomic_weight_ratio == pytest.approx(236.9986)


def test_fission_energy(pu239):
    fer = pu239.fission_energy
    assert isinstance(fer, openmc.data.FissionEnergyRelease)
    components = ['betas', 'delayed_neutrons', 'delayed_photons', 'fragments',
                  'neutrinos', 'prompt_neutrons', 'prompt_photons', 'recoverable',
                  'total', 'q_prompt', 'q_recoverable', 'q_total']
    for c in components:
        assert isinstance(getattr(fer, c), Callable)


def test_compact_fission_energy(tmpdir):
    files = [os.path.join(_ENDF_DATA, 'neutrons', 'n-090_Th_232.endf'),
             os.path.join(_ENDF_DATA, 'neutrons', 'n-094_Pu_240.endf'),
             os.path.join(_ENDF_DATA, 'neutrons', 'n-094_Pu_241.endf')]
    output = str(tmpdir.join('compact_lib.h5'))
    openmc.data.write_compact_458_library(files, output)
    assert os.path.exists(output)


def test_energy_grid(pu239):
    assert isinstance(pu239.energy, Mapping)
    for temp, grid in pu239.energy.items():
        assert temp.endswith('K')
        assert np.all(np.diff(grid) >= 0.0)


def test_reactions(pu239):
    assert 2 in pu239.reactions
    assert isinstance(pu239.reactions[2], openmc.data.Reaction)
    with pytest.raises(KeyError):
        pu239.reactions[14]


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


def test_derived_products(am244):
    fission = am244.reactions[18]
    total_neutron = fission.derived_products[0]
    assert total_neutron.emission_mode == 'total'
    assert total_neutron.yield_(6e6) == pytest.approx(4.2558)


def test_urr(pu239):
    for T, ptable in pu239.urr.items():
        assert T.endswith('K')
        assert isinstance(ptable, openmc.data.ProbabilityTables)
    ptable = pu239.urr['294K']
    assert ptable.absorption_flag == -1
    assert ptable.energy[0] == pytest.approx(2500.001)
    assert ptable.energy[-1] == pytest.approx(29999.99)
    assert ptable.inelastic_flag == 51
    assert ptable.interpolation == 2
    assert not ptable.multiply_smooth
    assert ptable.table.shape == (70, 6, 20)
    assert ptable.table.shape[0] == ptable.energy.size


def test_get_reaction_components(h2):
    assert h2.get_reaction_components(1) == [2, 16, 102]
    assert h2.get_reaction_components(101) == [102]
    assert h2.get_reaction_components(16) == [16]
    assert h2.get_reaction_components(51) == []


def test_export_to_hdf5(tmpdir, pu239, gd154):
    filename = str(tmpdir.join('pu239.h5'))
    pu239.export_to_hdf5(filename)
    assert os.path.exists(filename)
    with pytest.raises(NotImplementedError):
        gd154.export_to_hdf5('gd154.h5')

def test_slbw(xe135):
    res = xe135.resonances
    assert isinstance(res, openmc.data.Resonances)
    assert len(res.ranges) == 2
    resolved = res.resolved
    assert isinstance(resolved, openmc.data.SingleLevelBreitWigner)
    assert resolved.energy_min == pytest.approx(1e-5)
    assert resolved.energy_max == pytest.approx(190.)
    assert resolved.target_spin == pytest.approx(1.5)
    assert isinstance(resolved.parameters, pd.DataFrame)
    s = resolved.parameters.iloc[0]
    assert s['energy'] == pytest.approx(0.084)

    xs = resolved.reconstruct([10., 30., 100.])
    assert sorted(xs.keys()) == [2, 18, 102]
    assert np.all(xs[18] == 0.0)


def test_mlbw(sm150):
    resolved = sm150.resonances.resolved
    assert isinstance(resolved, openmc.data.MultiLevelBreitWigner)
    assert resolved.energy_min == pytest.approx(1e-5)
    assert resolved.energy_max == pytest.approx(1570.)
    assert resolved.target_spin == 0.0

    xs = resolved.reconstruct([10., 100., 1000.])
    assert sorted(xs.keys()) == [2, 18, 102]
    assert np.all(xs[18] == 0.0)


def test_reichmoore(gd154):
    res = gd154.resonances
    assert isinstance(res, openmc.data.Resonances)
    assert len(res.ranges) == 2
    resolved, unresolved = res.ranges
    assert resolved is res.resolved
    assert unresolved is res.unresolved
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

    elastic = gd154.reactions[2].xs['0K']
    assert isinstance(elastic, openmc.data.ResonancesWithBackground)
    assert elastic(0.0253) == pytest.approx(5.7228949796394524)


def test_rml(cl35):
    resolved = cl35.resonances.resolved
    assert isinstance(resolved, openmc.data.RMatrixLimited)
    assert resolved.energy_min == pytest.approx(1e-5)
    assert resolved.energy_max == pytest.approx(1.2e6)
    assert resolved.target_spin == 0.0
    for group in resolved.spin_groups:
        assert isinstance(group, openmc.data.SpinGroup)


def test_madland_nix(am241):
    fission = am241.reactions[18]
    prompt_neutron = fission.products[0]
    dist = prompt_neutron.distribution[0].energy
    assert isinstance(dist, openmc.data.MadlandNix)
    assert dist.efl == pytest.approx(1029979.0)
    assert dist.efh == pytest.approx(546729.7)
    assert isinstance(dist.tm, Callable)


def test_watt(u233):
    fission = u233.reactions[18]
    prompt_neutron = fission.products[0]
    dist = prompt_neutron.distribution[0].energy
    assert isinstance(dist, openmc.data.WattEnergy)


def test_maxwell(u236):
    fission = u236.reactions[18]
    prompt_neutron = fission.products[0]
    dist = prompt_neutron.distribution[0].energy
    assert isinstance(dist, openmc.data.MaxwellEnergy)


def test_evaporation(na22):
    n2n = na22.reactions[16]
    dist = n2n.products[0].distribution[0].energy
    assert isinstance(dist, openmc.data.Evaporation)


def test_laboratory(be9):
    n2n = be9.reactions[16]
    dist = n2n.products[0].distribution[0]
    assert isinstance(dist, openmc.data.LaboratoryAngleEnergy)
    assert list(dist.breakpoints) == [18]
    assert list(dist.interpolation) == [2]
    assert dist.energy[0] == pytest.approx(1748830.)
    assert dist.energy[-1] == pytest.approx(20.e6)
    assert len(dist.energy) == len(dist.energy_out) == len(dist.mu)
    for eout, mu in zip(dist.energy_out, dist.mu):
        assert len(eout) == len(mu)
        assert np.all((-1. <= mu.x) & (mu.x <= 1.))


def test_correlated(tmpdir):
    endf_file = os.path.join(_ENDF_DATA, 'neutrons', 'n-014_Si_030.endf')
    si30 = openmc.data.IncidentNeutron.from_njoy(endf_file, heatr=False)

    # Convert to HDF5 and read back
    filename = str(tmpdir.join('si30.h5'))
    si30.export_to_hdf5(filename)
    si30_copy = openmc.data.IncidentNeutron.from_hdf5(filename)


def test_nbody(tmpdir, h2):
    # Convert to HDF5 and read back
    filename = str(tmpdir.join('h2.h5'))
    h2.export_to_hdf5(filename)
    h2_copy = openmc.data.IncidentNeutron.from_hdf5(filename)

    # Compare distributions
    nbody1 = h2[16].products[0].distribution[0]
    nbody2 = h2_copy[16].products[0].distribution[0]
    assert nbody1.total_mass == nbody2.total_mass
    assert nbody1.n_particles == nbody2.n_particles
    assert nbody1.q_value == nbody2.q_value


def test_ace_convert(tmpdir):
    filename = os.path.join(_ENDF_DATA, 'neutrons', 'n-001_H_001.endf')
    ace_ascii = str(tmpdir.join('ace_ascii'))
    ace_binary = str(tmpdir.join('ace_binary'))
    openmc.data.njoy.make_ace(filename, ace=ace_ascii)

    # Convert to binary
    openmc.data.ace.ascii_to_binary(ace_ascii, ace_binary)

    # Make sure conversion worked
    lib_ascii = openmc.data.ace.Library(ace_ascii)
    lib_binary = openmc.data.ace.Library(ace_binary)
    for tab_a, tab_b in zip(lib_ascii.tables, lib_binary.tables):
        assert tab_a.name == tab_b.name
        assert tab_a.atomic_weight_ratio == pytest.approx(tab_b.atomic_weight_ratio)
        assert tab_a.temperature == pytest.approx(tab_b.temperature)
        assert np.all(tab_a.nxs == tab_b.nxs)
        assert np.all(tab_a.jxs == tab_b.jxs)
