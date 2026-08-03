from collections.abc import Mapping, Callable
import inspect
import os

import numpy as np
import pandas as pd
import pytest
import h5py
import openmc.data
import openmc.data.neutron as neutron_module

from . import needs_njoy

_TEMPERATURES = [300., 600., 900.]


@pytest.fixture(scope='module')
def pu239():
    """Pu239 HDF5 data."""
    directory = os.path.dirname(openmc.config.get('cross_sections'))
    filename = os.path.join(directory, 'Pu239.h5')
    return openmc.data.IncidentNeutron.from_hdf5(filename)


@pytest.fixture(scope='module')
def xe135(endf_data):
    """Xe135 ENDF data (contains SLBW resonance range)"""
    filename = os.path.join(endf_data, 'neutrons', 'n-054_Xe_135.endf')
    return openmc.data.IncidentNeutron.from_endf(filename)


@pytest.fixture(scope='module')
def sm150(endf_data):
    """Sm150 ENDF data (contains MLBW resonance range)"""
    filename = os.path.join(endf_data, 'neutrons', 'n-062_Sm_150.endf')
    return openmc.data.IncidentNeutron.from_endf(filename)


@pytest.fixture(scope='module')
def gd154(endf_data):
    """Gd154 ENDF data (contains Reich Moore resonance range and reosnance
    covariance with LCOMP=1)."""
    filename = os.path.join(endf_data, 'neutrons', 'n-064_Gd_154.endf')
    return openmc.data.IncidentNeutron.from_endf(filename, covariance=True)


@pytest.fixture(scope='module')
def cl35(endf_data):
    """Cl35 ENDF data (contains RML resonance range)"""
    filename = os.path.join(endf_data, 'neutrons', 'n-017_Cl_035.endf')
    return openmc.data.IncidentNeutron.from_endf(filename)


@pytest.fixture(scope='module')
def am241(endf_data):
    """Am241 ENDF data (contains Madland-Nix fission energy distribution)."""
    filename = os.path.join(endf_data, 'neutrons', 'n-095_Am_241.endf')
    return openmc.data.IncidentNeutron.from_endf(filename)


@pytest.fixture(scope='module')
def u233(endf_data):
    """U233 ENDF data (contains Watt fission energy distribution)."""
    filename = os.path.join(endf_data, 'neutrons', 'n-092_U_233.endf')
    return openmc.data.IncidentNeutron.from_endf(filename)


@pytest.fixture(scope='module')
def u236(endf_data):
    """U236 ENDF data (contains Watt fission energy distribution)."""
    filename = os.path.join(endf_data, 'neutrons', 'n-092_U_236.endf')
    return openmc.data.IncidentNeutron.from_endf(filename)


@pytest.fixture(scope='module')
def na22(endf_data):
    """Na22 ENDF data (contains evaporation spectrum)."""
    filename = os.path.join(endf_data, 'neutrons', 'n-011_Na_022.endf')
    return openmc.data.IncidentNeutron.from_endf(filename)


@pytest.fixture(scope='module')
def na23(endf_data):
    """Na23 ENDF data (contains MLBW resonance covariance with LCOMP=0)."""
    filename = os.path.join(endf_data, 'neutrons', 'n-011_Na_023.endf')
    return openmc.data.IncidentNeutron.from_endf(filename, covariance=True)


@pytest.fixture(scope='module')
def be9(endf_data):
    """Be9 ENDF data (contains laboratory angle-energy distribution)."""
    filename = os.path.join(endf_data, 'neutrons', 'n-004_Be_009.endf')
    return openmc.data.IncidentNeutron.from_endf(filename)


@pytest.fixture(scope='module')
def h2(endf_data):
    endf_file = os.path.join(endf_data, 'neutrons', 'n-001_H_002.endf')
    return openmc.data.IncidentNeutron.from_njoy(
        endf_file, temperatures=_TEMPERATURES)


@pytest.fixture(scope='module')
def am244(endf_data):
    endf_file = os.path.join(endf_data, 'neutrons', 'n-095_Am_244.endf')
    return openmc.data.IncidentNeutron.from_njoy(endf_file)


@pytest.fixture(scope='module')
def ti50(endf_data):
    """Ti50 ENDF data (contains Multi-level Breit-Wigner resonance range and
       resonance covariance with LCOMP=1)."""
    filename = os.path.join(endf_data, 'neutrons', 'n-022_Ti_050.endf')
    return openmc.data.IncidentNeutron.from_endf(filename, covariance=True)


@pytest.fixture(scope='module')
def cf252(endf_data):
    """Cf252 ENDF data (contains RM resonance covariance with LCOMP=0)."""
    filename = os.path.join(endf_data, 'neutrons', 'n-098_Cf_252.endf')
    return openmc.data.IncidentNeutron.from_endf(filename, covariance=True)


@pytest.fixture(scope='module')
def th232(endf_data):
    """Th232 ENDF data (contains RM resonance covariance with LCOMP=2)."""
    filename = os.path.join(endf_data, 'neutrons', 'n-090_Th_232.endf')
    return openmc.data.IncidentNeutron.from_endf(filename, covariance=True)


def test_attributes(pu239):
    assert pu239.name == 'Pu239'
    assert pu239.mass_number == 239
    assert pu239.metastable == 0
    assert pu239.atomic_symbol == 'Pu'
    assert pu239.atomic_weight_ratio == pytest.approx(236.9986)


def test_from_endf_material(endf_data):
    filename = os.path.join(endf_data, 'neutrons', 'n-001_H_001.endf')
    material = openmc.data.endf.get_evaluations(filename)[0]

    data = openmc.data.IncidentNeutron.from_endf(material)

    assert data.name == 'H1'
    assert data.atomic_number == 1
    assert data.mass_number == 1
    assert 2 in data.reactions


def test_fission_energy_from_endf_material(endf_data):
    filename = os.path.join(endf_data, 'neutrons', 'n-092_U_235.endf')
    material = openmc.data.endf.get_evaluations(filename)[0]
    neutron = openmc.data.IncidentNeutron.from_endf(material)

    data = openmc.data.FissionEnergyRelease.from_endf(material, neutron)

    assert data.fragments(0.0) > 0.0


def test_fission_energy(pu239):
    fer = pu239.fission_energy
    assert isinstance(fer, openmc.data.FissionEnergyRelease)
    components = ['betas', 'delayed_neutrons', 'delayed_photons', 'fragments',
                  'neutrinos', 'prompt_neutrons', 'prompt_photons', 'recoverable',
                  'total', 'q_prompt', 'q_recoverable', 'q_total']
    for c in components:
        assert isinstance(getattr(fer, c), Callable)


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


@needs_njoy
def test_derived_products(am244):
    fission = am244.reactions[18]
    total_neutron = fission.derived_products[0]
    assert total_neutron.emission_mode == 'total'
    assert total_neutron.yield_(6e6) == pytest.approx(4.2558)


@needs_njoy
def test_kerma(run_in_tmpdir, am244, h2):
    # Make sure kerma w/ local photon is >= regular kerma
    for nuc in (am244, h2):
        assert 301 in nuc
        assert 901 in nuc
        for T in nuc.temperatures:
            k, k_local = nuc[301].xs[T], nuc[901].xs[T]
            assert np.all(k.x == k_local.x)
            assert np.all(k_local.y >= k.y)

    # Make sure 301/901 get exported/imported correctly
    h2.export_to_hdf5("H2.h5")
    read_in = openmc.data.IncidentNeutron.from_hdf5("H2.h5")
    assert 301 in read_in
    assert 901 in read_in
    assert np.all(read_in[901].xs['300K'].y == h2[901].xs['300K'].y)


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


@needs_njoy
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


def test_mlbw(sm150):
    resolved = sm150.resonances.resolved
    assert isinstance(resolved, openmc.data.MultiLevelBreitWigner)
    assert resolved.energy_min == pytest.approx(1e-5)
    assert resolved.energy_max == pytest.approx(1570.)
    assert resolved.target_spin == 0.0


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


def test_rml(cl35):
    resolved = cl35.resonances.resolved
    assert isinstance(resolved, openmc.data.RMatrixLimited)
    assert resolved.energy_min == pytest.approx(1e-5)
    assert resolved.energy_max == pytest.approx(1.2e6)
    assert resolved.target_spin == 0.0
    for group in resolved.spin_groups:
        assert isinstance(group, openmc.data.SpinGroup)


def test_mlbw_cov_lcomp0(cf252):
    # Testing on first range only
    cov = cf252.resonance_covariance.ranges[0]
    res = cf252.resonances.ranges[0]
    assert cov.parameters['energy'][0] == pytest.approx(-3.5)
    assert res.parameters['energy'][0] == cov.parameters['energy'][0]
    assert isinstance(cov, openmc.data.resonance_covariance.MultiLevelBreitWignerCovariance)
    assert cov.energy_min == pytest.approx(1e-5)
    assert cov.energy_max == pytest.approx(1000.)
    assert cov.covariance[0,0] == pytest.approx(1.225e-05)

    subset = cov.subset('energy', [0, 100])
    assert not subset.parameters.empty
    assert (subset.file2res.parameters['energy'] < 100).all()
    samples = cov.sample(1)


def test_mlbw_cov_lcomp1(ti50):
    # Testing on first range only
    cov = ti50.resonance_covariance.ranges[0]
    res = ti50.resonances.ranges[0]
    assert cov.parameters['energy'][0] == pytest.approx(-21020.)
    assert res.parameters['energy'][0] == cov.parameters['energy'][0]
    assert isinstance(cov, openmc.data.resonance_covariance.MultiLevelBreitWignerCovariance)
    assert cov.energy_min == pytest.approx(1e-5)
    assert cov.energy_max == pytest.approx(587000.)
    assert cov.covariance[0,0] == pytest.approx(1.410177e5)

    subset = cov.subset('L', [1, 1])
    assert not subset.parameters.empty
    assert (subset.file2res.parameters['L'] == 1).all()
    cov.sample(1)


def test_mlbw_cov_lcomp2(na23):
    # Testing on first range only
    cov = na23.resonance_covariance.ranges[0]
    res = na23.resonances.ranges[0]
    assert cov.parameters['energy'][0] == pytest.approx(2810.)
    assert res.parameters['energy'][0] == cov.parameters['energy'][0]
    assert isinstance(cov, openmc.data.resonance_covariance.MultiLevelBreitWignerCovariance)
    assert cov.energy_min == pytest.approx(600)
    assert cov.energy_max == pytest.approx(500000.)
    assert cov.covariance[0,0] == pytest.approx(16.1064163584)

    subset = cov.subset('L', [1, 1])
    assert not subset.parameters.empty
    assert (subset.file2res.parameters['L'] == 1).all()
    cov.sample(1)


def test_rmcov_lcomp1(gd154):
    # Testing on first range only
    cov = gd154.resonance_covariance.ranges[0]
    res = gd154.resonances.ranges[0]
    assert cov.parameters['energy'][0] == pytest.approx(-2.200001)
    assert res.parameters['energy'][0] == cov.parameters['energy'][0]
    assert isinstance(cov, openmc.data.resonance_covariance.ReichMooreCovariance)
    assert cov.energy_min == pytest.approx(1e-5)
    assert cov.energy_max == pytest.approx(2760.)
    assert cov.covariance[0,0] == pytest.approx(0.8895997)

    subset = cov.subset('energy', [0, 100])
    assert not subset.parameters.empty
    assert (subset.file2res.parameters['energy'] < 100).all()
    cov.sample(1)


def test_rmcov_lcomp2(th232):
    # Testing on first range only
    cov = th232.resonance_covariance.ranges[0]
    res = th232.resonances.ranges[0]
    assert cov.parameters['energy'][0] == pytest.approx(-2000)
    assert res.parameters['energy'][0] == cov.parameters['energy'][0]
    assert isinstance(cov, openmc.data.resonance_covariance.ReichMooreCovariance)
    assert cov.energy_min == pytest.approx(1e-5)
    assert cov.energy_max == pytest.approx(4000.)
    assert cov.covariance[0,0] == pytest.approx(246.6043092496)

    subset = cov.subset('energy', [0, 100])
    assert not subset.parameters.empty
    assert (subset.file2res.parameters['energy'] < 100).all()
    cov.sample(1)


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


@needs_njoy
def test_correlated(tmpdir, endf_data):
    endf_file = os.path.join(endf_data, 'neutrons', 'n-014_Si_030.endf')
    si30 = openmc.data.IncidentNeutron.from_njoy(endf_file, heatr=False)

    # Convert to HDF5 and read back
    filename = str(tmpdir.join('si30.h5'))
    si30.export_to_hdf5(filename)
    si30_copy = openmc.data.IncidentNeutron.from_hdf5(filename)


@needs_njoy
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


@needs_njoy
def test_ace_convert(run_in_tmpdir, endf_data):
    filename = os.path.join(endf_data, 'neutrons', 'n-001_H_001.endf')
    ace_ascii = 'ace_ascii'
    ace_binary = 'ace_binary'
    openmc.data.njoy.make_ace(filename, acer=ace_ascii)

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
        assert tab_a.zaid == tab_b.zaid
        assert tab_a.data_type == tab_b.data_type


def test_ace_table_types():
    TT = openmc.data.ace.TableType
    assert TT.from_suffix('c') == TT.NEUTRON_CONTINUOUS
    assert TT.from_suffix('nc') == TT.NEUTRON_CONTINUOUS
    assert TT.from_suffix('80c') == TT.NEUTRON_CONTINUOUS
    assert TT.from_suffix('t') == TT.THERMAL_SCATTERING
    assert TT.from_suffix('20t') == TT.THERMAL_SCATTERING
    with pytest.raises(ValueError):
        TT.from_suffix('z')


@needs_njoy
def test_high_temperature(endf_data):
    endf_file = os.path.join(endf_data, 'neutrons', 'n-001_H_001.endf')

    # Ensure that from_njoy works when given a high temperature
    openmc.data.IncidentNeutron.from_njoy(endf_file, temperatures=[123_456.0])


class _MetadataEvaluation:
    def __init__(self, text, mt=51):
        self.section = {(6, mt): text}


def _endf_record(c1=0.0, c2=0.0, l1=0, l2=0, n1=0, n2=0):
    """Minimal fixed-width ENDF record for direct MF=6 parser tests."""
    return (f'{c1:11.5E}{c2:11.5E}{l1:11d}{l2:11d}'
            f'{n1:11d}{n2:11d}\n')


def _endf_values(*values):
    """One or more 66-column ENDF floating-point payload records."""
    records = []
    for start in range(0, len(values), 6):
        fields = ''.join(f'{value:11.5E}' for value in values[start:start + 6])
        records.append(f'{fields:66s}\n')
    return ''.join(records)


def _endf_integers(*values):
    fields = ''.join(f'{value:11d}' for value in values)
    return f'{fields:66s}\n'


def _metadata_tab1(zap, awp, lip, law):
    return (_endf_record(zap, awp, lip, law, 1, 2) +
            _endf_integers(2, 2) +
            _endf_values(0.0, 0.0, 20.0e6, 1.0))


def _metadata_tab2(l1=0, n2=0):
    return (_endf_record(0.0, 0.0, l1, 0, 1, n2) +
            _endf_integers(n2, 2))


def _metadata_list(l1, l2, n1, n2):
    return (_endf_record(0.0, 0.0, l1, l2, n1, n2) +
            _endf_values(*range(n1)))


def _metadata_section(mt, lct, entries, jp=0, trailing=''):
    text = _endf_record(26056.0, 55.0, jp, lct, len(entries), 0)
    for zap, awp, lip, law, payload in entries:
        text += _metadata_tab1(zap, awp, lip, law) + payload
    return _MetadataEvaluation(text + trailing, mt)


def _law1(lang, lists):
    payload = _metadata_tab2(lang, len(lists))
    for nep, na, nd, nw in lists:
        payload += _metadata_list(nd, na, nw, nep)
    return payload


def _law2(lists):
    payload = _metadata_tab2(0, len(lists))
    for lang, nl, nw in lists:
        payload += _metadata_list(lang, 0, nw, nl)
    return payload


def _law5(lidp, lists):
    payload = _metadata_tab2(lidp, len(lists))
    for ltp, nl, nw in lists:
        payload += _metadata_list(ltp, 0, nw, nl)
    return payload


def _law6(apsx=4.0, c2=0.0, l1=0, l2=0, n1=0, npsx=3):
    npsx_field = (f'{npsx:11d}' if isinstance(npsx, int)
                  else f'{npsx:>11}')
    return (f'{apsx:11.5E}{c2:11.5E}{l1:11d}{l2:11d}'
            f'{n1:11d}{npsx_field}\n')


def _law7(nmu_values):
    payload = _metadata_tab2(0, len(nmu_values))
    for nmu in nmu_values:
        payload += _metadata_tab2(0, nmu)
        payload += _metadata_tab1(0, 0.0, 0, 0) * nmu
    return payload


def _matrix_cases():
    law1 = (1, 1.0, 0, 1,
            _law1(2, [(3, 1, 0, 9), (2, 2, 2, 8)]))
    law1_bad = (1, 1.0, 0, 1, _law1(1, [(2, 0, 0, 3)]))
    law1_even = (1, 1.0, 0, 1, _law1(12, [(2, 2, 1, 8)]))
    law2 = (1, 1.0, 0, 2, _law2([(0, 3, 3), (14, 2, 4)]))
    capture = (0, 6129000.0, 0, 2, _law2([(0, 1, 1)]))
    law5 = (1001, 1.0, 0, 5,
            _law5(1, [(1, 2, 9), (2, 2, 3), (15, 2, 4)]))
    law6 = (1, 1.0, 0, 6, _law6())
    law7 = (1, 1.0, 0, 7, _law7([0, 2]))
    fission_negative = [
        (1, 1.0, 0, 0, ''), (1, 1.0, 0, -5, ''),
        (0, 0.0, 0, 0, ''), (0, 0.0, 0, -15, '')]
    fission_positive = [
        (1, 1.0, 0, 0, ''), (1, 1.0, 0, 3, ''),
        (0, 0.0, 0, 0, ''), (0, 0.0, 0, 4, '')]
    return [
        ('law1_exact_counts', _metadata_section(51, 3, [law1]),
         'published_structurally_complete'),
        ('law1_bad_nw', _metadata_section(51, 1, [law1_bad]),
         'count_mismatch'),
        ('law1_tabulated_even_na', _metadata_section(51, 1, [law1_even]),
         'published_structurally_complete'),
        ('law2_cm_and_counts', _metadata_section(51, 2, [law2]),
         'published_structurally_complete'),
        ('law2_capture_photon_awp_ev', _metadata_section(102, 2, [capture]),
         'published_structurally_complete'),
        ('law5_exact_lidp_counts', _metadata_section(51, 1, [law5]),
         'published_structurally_complete'),
        ('law6_controls', _metadata_section(51, 1, [law6]),
         'published_structurally_complete'),
        ('law7_nested_counts', _metadata_section(51, 1, [law7]),
         'published_structurally_complete'),
        ('law0_zero_payload', _metadata_section(
            51, 1, [(1, 1.0, 0, 0, '')]),
         'published_structurally_complete'),
        ('law3_zero_payload', _metadata_section(
            51, 1, [(1, 1.0, 0, 3, '')]),
         'published_structurally_complete'),
        ('law4_zero_payload', _metadata_section(
            51, 1, [(1, 1.0, 0, 4, '')]),
         'published_structurally_complete'),
        ('negative_law_fission_order', _metadata_section(
            18, 1, fission_negative, jp=11),
         'published_structurally_complete'),
        ('jp2_positive_law_sequence', _metadata_section(
            18, 1, fission_positive, jp=22),
         'published_structurally_complete'),
        ('unknown_law_never_publishes', _metadata_section(
            51, 1, [(1, 1.0, 0, 99, '')]), 'unknown_law'),
        ('trailing_record_never_publishes', _metadata_section(
            51, 1, [(1, 1.0, 0, 0, '')], trailing='unexpected'),
         'trailing_record'),
    ]


@pytest.mark.parametrize('name,evaluation,expected', _matrix_cases(),
                         ids=[case[0] for case in _matrix_cases()])
def test_evaluated_product_metadata_matrix(name, evaluation, expected):
    mt = next(mt for mf, mt in evaluation.section if mf == 6)
    result = openmc.data.get_evaluated_product_metadata(evaluation, mt)
    if expected == 'published_structurally_complete':
        assert result.result_status == expected
        assert result.structural_failure is None
        assert result.metadata is not None
        assert result.hdf5_group_version == 1
        assert result.entries == result.metadata.entries
        assert (result.jp, result.lct, result.section_semantic_status) == (
            result.metadata.jp, result.metadata.lct,
            result.metadata.section_semantic_status)
        if name == 'law2_capture_photon_awp_ev':
            entry = result.metadata.entries[0]
            assert entry.awp == 6129000.0
            assert entry.awp_interpretation == 'primary_photon_energy_eV'
        elif name == 'law5_exact_lidp_counts':
            entry = result.metadata.entries[0]
            assert (entry.identity_status, entry.derived_particle) == (
                'raw_only', None)
            assert entry.ltp_values == (1, 2, 15)
        elif name == 'law6_controls':
            entry = result.metadata.entries[0]
            assert (entry.apsx, entry.npsx) == (4.0, 3)
        elif name == 'negative_law_fission_order':
            assert tuple(entry.law for entry in result.metadata.entries) == (
                0, -5, 0, -15)
    else:
        assert result.result_status == 'structural_failure'
        assert result.structural_failure == expected
        assert result.metadata is None
        assert (result.jp, result.lct, result.section_semantic_status,
                result.hdf5_group_version, result.entries) == (
            None, None, None, None, ())


def test_evaluated_product_metadata_direct_parser():
    law0 = _metadata_section(51, 1, [(1, 1.0, 0, 0, '')])
    result = openmc.data.get_evaluated_product_metadata(law0, 51)
    assert result.result_status == 'published_structurally_complete'
    assert result.metadata.entries[0].subsection_index == 0
    assert result.metadata.entries[0].derived_particle == 'neutron'
    assert result.metadata.entries[0].law == 0
    assert (result.metadata.entries[0].apsx,
            result.metadata.entries[0].npsx) == (None, None)

    with pytest.raises(ValueError):
        openmc.data.EvaluatedProductMetadataEntry(
            0, 1, 1.0, 'mass_ratio', 0, 6, (), None, (), 'complete',
            'special_particle', 'neutron')
    for invalid_apsx in (True, '4.0'):
        with pytest.raises(ValueError):
            openmc.data.EvaluatedProductMetadataEntry(
                0, 1, 1.0, 'mass_ratio', 0, 6, (), None, (), 'complete',
                'special_particle', 'neutron', invalid_apsx, 3)
    with pytest.raises(ValueError):
        openmc.data.EvaluatedProductMetadataEntry(
            0, 1, 1.0, 'mass_ratio', 0, 0, (), None, (), 'complete',
            'special_particle', 'neutron', 4.0, 3)

    absent = openmc.data.get_evaluated_product_metadata(
        _MetadataEvaluation('', 52), 51)
    assert absent == openmc.data.EvaluatedProductMetadataResult(51, 'absent')
    assert (absent.schema, absent.source_format, absent.mf, absent.mt) == (
        'openmc-evaluated-product-metadata/v2', 'ENDF-6', 6, 51)
    assert (absent.jp, absent.lct, absent.section_semantic_status,
            absent.hdf5_group_version, absent.entries) == (
        None, None, None, None, ())

    premature = _MetadataEvaluation(
        _endf_record(26056.0, 55.0, 0, 1, 1, 0))
    result = openmc.data.get_evaluated_product_metadata(premature, 51)
    assert result.structural_failure == 'premature_end'

    empty_active = _MetadataEvaluation(
        _endf_record(26056.0, 55.0, 1, 1, 0, 0), 18)
    result = openmc.data.get_evaluated_product_metadata(empty_active, 18)
    assert result.result_status == 'structural_failure'
    assert result.structural_failure == 'invalid_header'

    nonzero_head_n2 = _MetadataEvaluation(
        _endf_record(26056.0, 55.0, 0, 1, 0, 1))
    result = openmc.data.get_evaluated_product_metadata(nonzero_head_n2, 51)
    assert result.structural_failure == 'invalid_header'


def test_evaluated_product_metadata_semantic_refusals():
    invalid_lidp = _metadata_section(51, 1, [
        (1001, 1.0, 0, 5, _law5(2, [(1, 2, 1)]))])
    result = openmc.data.get_evaluated_product_metadata(invalid_lidp, 51)
    entry = result.metadata.entries[0]
    assert entry.identity_status == 'raw_only'
    assert entry.derived_particle is None
    assert entry.semantic_status == 'invalid_combination'

    unsupported = _metadata_section(51, 1, [
        (1, 1.0, 0, 1, _law1(999, [(1, 0, 0, 2)]))])
    result = openmc.data.get_evaluated_product_metadata(unsupported, 51)
    assert result.metadata.section_semantic_status == 'unsupported_lang'

    inactive_negative = _metadata_section(
        18, 1, [(1, 1.0, 0, -5, '')], jp=0)
    result = openmc.data.get_evaluated_product_metadata(inactive_negative, 18)
    assert result.metadata.section_semantic_status == 'invalid_combination'

    chance_fission = _metadata_section(19, 1, [
        (1, 1.0, 0, 0, ''), (1, 1.0, 0, -5, '')], jp=1)
    result = openmc.data.get_evaluated_product_metadata(chance_fission, 19)
    assert result.metadata.section_semantic_status == 'invalid_combination'

    crossed_neutron = _metadata_section(18, 1, [
        (1, 1.0, 0, -15, ''),
        (0, 0.0, 0, 0, ''), (0, 0.0, 0, -15, '')], jp=10)
    result = openmc.data.get_evaluated_product_metadata(crossed_neutron, 18)
    assert result.metadata.section_semantic_status == 'invalid_combination'

    crossed_photon = _metadata_section(18, 1, [
        (1, 1.0, 0, 0, ''), (1, 1.0, 0, -5, ''),
        (0, 0.0, 0, -5, '')], jp=1)
    result = openmc.data.get_evaluated_product_metadata(crossed_photon, 18)
    assert result.metadata.section_semantic_status == 'invalid_combination'

    bad_tabulated_na = _metadata_section(51, 1, [
        (1, 1.0, 0, 1, _law1(12, [(1, 0, 0, 2)]))])
    result = openmc.data.get_evaluated_product_metadata(bad_tabulated_na, 51)
    assert result.structural_failure == 'count_mismatch'


@pytest.mark.parametrize('lang,na,expected', [
    (1, 0, 'complete'),
    (2, 0, 'complete'),
    (2, 1, 'complete'),
    (2, 2, 'complete'),
    (11, 2, 'complete'),
    (12, 2, 'complete'),
    (13, 2, 'complete'),
    (14, 2, 'complete'),
    (15, 2, 'complete'),
    (999, 0, 'unsupported_lang'),
])
def test_evaluated_product_metadata_law1_selectors(lang, na, expected):
    nw = na + 2
    evaluation = _metadata_section(51, 2 if lang == 2 else 1, [
        (1, 1.0, 0, 1, _law1(lang, [(1, na, 0, nw)]))])
    result = openmc.data.get_evaluated_product_metadata(evaluation, 51)
    assert result.result_status == 'published_structurally_complete'
    assert result.entries[0].semantic_status == expected


@pytest.mark.parametrize('lct,expected', [
    (1, 'invalid_combination'),
    (2, 'complete'),
    (3, 'invalid_combination'),
    (4, 'invalid_combination'),
])
def test_evaluated_product_metadata_law1_lang2_lct_semantics(lct, expected):
    evaluation = _metadata_section(51, lct, [
        (1, 1.0, 0, 1, _law1(2, [(1, 0, 0, 2)]))])
    result = openmc.data.get_evaluated_product_metadata(evaluation, 51)

    assert result.result_status == 'published_structurally_complete'
    assert result.metadata is not None
    assert result.entries[0].semantic_status == expected
    assert result.section_semantic_status == expected


@pytest.mark.parametrize('lang,nl,nw,expected', [
    (0, 1, 1, 'complete'),
    (12, 2, 4, 'complete'),
    (14, 2, 4, 'complete'),
    (999, 0, 0, 'unsupported_lang'),
])
def test_evaluated_product_metadata_law2_selectors(lang, nl, nw, expected):
    evaluation = _metadata_section(51, 2, [
        (1, 1.0, 0, 2, _law2([(lang, nl, nw)]))])
    result = openmc.data.get_evaluated_product_metadata(evaluation, 51)
    assert result.result_status == 'published_structurally_complete'
    assert result.entries[0].semantic_status == expected


@pytest.mark.parametrize('lidp,ltp,nl,nw,expected', [
    (0, 1, 1, 7, 'complete'),
    (1, 1, 1, 6, 'complete'),
    (0, 2, 1, 2, 'complete'),
    (0, 12, 1, 2, 'complete'),
    (0, 14, 1, 2, 'complete'),
    (0, 15, 1, 2, 'complete'),
    (0, 99, 0, 0, 'unsupported_ltp'),
    (2, 1, 1, 1, 'invalid_combination'),
])
def test_evaluated_product_metadata_law5_selectors(
        lidp, ltp, nl, nw, expected):
    evaluation = _metadata_section(51, 1, [
        (1001, 1.0, 0, 5, _law5(lidp, [(ltp, nl, nw)]))])
    result = openmc.data.get_evaluated_product_metadata(evaluation, 51)
    assert result.result_status == 'published_structurally_complete'
    assert result.entries[0].semantic_status == expected


def test_evaluated_product_metadata_law_envelope_mutations():
    missing_law6 = _metadata_section(51, 1, [(1, 1.0, 0, 6, '')])
    result = openmc.data.get_evaluated_product_metadata(missing_law6, 51)
    assert result.structural_failure == 'truncated_record'

    extra_law6 = _metadata_section(51, 1, [
        (1, 1.0, 0, 6, _law6())],
        trailing=_endf_record())
    result = openmc.data.get_evaluated_product_metadata(extra_law6, 51)
    assert result.structural_failure == 'trailing_record'

    for nmu_values in ([], [1], [0, 1, 2]):
        law7 = _metadata_section(51, 1, [
            (1, 1.0, 0, 7, _law7(nmu_values))])
        result = openmc.data.get_evaluated_product_metadata(law7, 51)
        assert result.result_status == 'published_structurally_complete'

    truncated_law7 = _metadata_section(51, 1, [
        (1, 1.0, 0, 7, _metadata_tab2(0, 1) + _metadata_tab2(0, 1))])
    result = openmc.data.get_evaluated_product_metadata(truncated_law7, 51)
    assert result.structural_failure == 'truncated_record'

    bad_law5 = _metadata_section(51, 1, [
        (1001, 1.0, 0, 5, _law5(0, [(12, 1, 1)]))])
    result = openmc.data.get_evaluated_product_metadata(bad_law5, 51)
    assert result.structural_failure == 'count_mismatch'


@pytest.mark.parametrize('payload', [
    _law6(0.0),
    _law6(float('nan')),
    _law6(4.0, c2=1.0),
    _law6(4.0, l1=1),
    _law6(4.0, l2=1),
    _law6(4.0, n1=1),
    _law6(4.0, npsx='3.5'),
    _law6(4.0, npsx=2),
    _law6(4.0, npsx=2**31),
], ids=[
    'nonpositive_apsx', 'nonfinite_apsx', 'reserved_c2', 'reserved_l1',
    'reserved_l2', 'reserved_n1', 'nonintegral_npsx', 'small_npsx',
    'unrepresentable_npsx',
])
def test_evaluated_product_metadata_law6_invalid_controls_are_atomic(payload):
    evaluation = _metadata_section(51, 1, [(1, 1.0, 0, 6, payload)])
    result = openmc.data.get_evaluated_product_metadata(evaluation, 51)

    assert result.result_status == 'structural_failure'
    assert result.structural_failure == 'invalid_law6_control'
    assert result.metadata is None
    assert result.entries == ()


def test_evaluated_product_metadata_law6_truncated_control_is_atomic():
    evaluation = _metadata_section(51, 1, [
        (1, 1.0, 0, 6, _law6()[:65])])
    result = openmc.data.get_evaluated_product_metadata(evaluation, 51)

    assert result.result_status == 'structural_failure'
    assert result.structural_failure == 'truncated_record'
    assert result.metadata is None
    assert result.entries == ()


def test_evaluated_product_metadata_boundaries_and_atomicity():
    lct4 = _metadata_section(51, 4, [(1, 1.0, 0, 0, '')])
    result = openmc.data.get_evaluated_product_metadata(lct4, 51)
    assert result.section_semantic_status == 'complete'

    lct3_nonlaw1 = _metadata_section(51, 3, [(1, 1.0, 0, 0, '')])
    result = openmc.data.get_evaluated_product_metadata(lct3_nonlaw1, 51)
    assert result.section_semantic_status == 'invalid_combination'

    unknown_lct = _MetadataEvaluation(
        _endf_record(26056.0, 55.0, 0, 5, 0, 0))
    result = openmc.data.get_evaluated_product_metadata(unknown_lct, 51)
    assert result.structural_failure == 'invalid_header'

    truncated_tab1 = _MetadataEvaluation(
        _endf_record(26056.0, 55.0, 0, 1, 1, 0) +
        _endf_record(1.0, 1.0, 0, 0, 1, 2))
    result = openmc.data.get_evaluated_product_metadata(truncated_tab1, 51)
    assert result.structural_failure == 'truncated_record'
    assert result.entries == ()

    truncated_list = _metadata_section(51, 1, [
        (1, 1.0, 0, 1,
         _metadata_tab2(1, 1) + _endf_record(0.0, 0.0, 0, 0, 2, 1))])
    result = openmc.data.get_evaluated_product_metadata(truncated_list, 51)
    assert result.structural_failure == 'truncated_record'
    assert result.entries == ()

    corrupt_list_length = _metadata_section(51, 1, [
        (1, 1.0, 0, 1,
         _metadata_tab2(1, 1) + _endf_record(0.0, 0.0, 0, 0, -1, 1))])
    result = openmc.data.get_evaluated_product_metadata(
        corrupt_list_length, 51)
    assert result.structural_failure == 'count_mismatch'
    assert result.entries == ()

    nk_over = _MetadataEvaluation(
        _endf_record(26056.0, 55.0, 0, 1, 2, 0) +
        _metadata_tab1(1, 1.0, 0, 0))
    result = openmc.data.get_evaluated_product_metadata(nk_over, 51)
    assert result.structural_failure == 'premature_end'

    nk_under = _MetadataEvaluation(
        _endf_record(26056.0, 55.0, 0, 1, 0, 0) +
        _metadata_tab1(1, 1.0, 0, 0))
    result = openmc.data.get_evaluated_product_metadata(nk_under, 51)
    assert result.structural_failure == 'trailing_record'


def test_evaluated_product_metadata_hdf5(tmp_path):
    entry = openmc.data.EvaluatedProductMetadataEntry(
        0, 0, 6129000.0, 'primary_photon_energy_eV', 0, 2, (0,), None,
        (), 'complete', 'special_particle', 'photon')
    metadata = openmc.data.EvaluatedProductMetadata(
        102, 0, 2, 'complete', (entry,))
    reaction = openmc.data.Reaction(102)
    transport_product = openmc.data.Product('neutron')
    reaction.products = [transport_product]
    reaction.evaluated_product_metadata = metadata
    filename = tmp_path / 'metadata.h5'
    with h5py.File(filename, 'w') as h5file:
        group = h5file.create_group('reaction_102')
        reaction.to_hdf5(group)
        metadata_group = group['evaluated_product_metadata']
        assert metadata_group.attrs['version'] == 1
        assert metadata_group.attrs['schema'] == (
            'openmc-evaluated-product-metadata/v2')
        assert metadata_group.attrs['result_status'] == (
            'published_structurally_complete')
        assert h5py.check_string_dtype(
            metadata_group.attrs.get_id('schema').dtype).length is None
        assert metadata_group.attrs.get_id('version').dtype == np.dtype('int32')
        assert set(metadata_group) == {'entry_0'}
        assert not any(name.startswith('product_') for name in metadata_group)
        assert metadata_group['entry_0']['lang_values'].dtype == np.dtype('int32')
        assert 'apsx' not in metadata_group['entry_0'].attrs
        assert 'npsx' not in metadata_group['entry_0'].attrs
        copied = openmc.data.Reaction.from_hdf5(group, {})
        assert copied.products[0].particle == 'neutron'
        assert copied.evaluated_product_metadata == metadata
        del group['evaluated_product_metadata']
        copied = openmc.data.Reaction.from_hdf5(group, {})
        assert copied.evaluated_product_metadata is None

    with h5py.File(filename, 'a') as h5file:
        reaction.to_hdf5(h5file.create_group('malformed'))
        h5file['malformed/evaluated_product_metadata'].attrs['version'] = 2
        with pytest.raises(ValueError):
            openmc.data.Reaction.from_hdf5(h5file['malformed'], {})

        fixed = h5file.create_group('fixed_string')
        reaction.to_hdf5(fixed)
        wire = fixed['evaluated_product_metadata']
        del wire.attrs['schema']
        wire.attrs['schema'] = np.bytes_(
            'openmc-evaluated-product-metadata/v2')
        with pytest.raises(ValueError):
            openmc.data.Reaction.from_hdf5(fixed, {})

        unknown_field = h5file.create_group('unknown_field')
        reaction.to_hdf5(unknown_field)
        unknown_field['evaluated_product_metadata'].attrs['extra'] = 1
        with pytest.raises(ValueError):
            openmc.data.Reaction.from_hdf5(unknown_field, {})

        gapped = h5file.create_group('gapped')
        reaction.to_hdf5(gapped)
        wire = gapped['evaluated_product_metadata']
        wire.move('entry_0', 'entry_1')
        with pytest.raises(ValueError):
            openmc.data.Reaction.from_hdf5(gapped, {})

        bad_rank = h5file.create_group('bad_rank')
        reaction.to_hdf5(bad_rank)
        entry_group = bad_rank['evaluated_product_metadata/entry_0']
        del entry_group['lang_values']
        entry_group.create_dataset(
            'lang_values', data=np.zeros((1, 1), dtype=np.int32))
        with pytest.raises(ValueError):
            openmc.data.Reaction.from_hdf5(bad_rank, {})

        inconsistent = h5file.create_group('inconsistent')
        reaction.to_hdf5(inconsistent)
        inconsistent['evaluated_product_metadata'].attrs['jp'] = np.int32(1)
        with pytest.raises(ValueError):
            openmc.data.Reaction.from_hdf5(inconsistent, {})

        portable = h5file.create_group('opposite_endian')
        reaction.to_hdf5(portable)
        wire = portable['evaluated_product_metadata']
        del wire.attrs['version']
        wire.attrs.create('version', np.asarray(1, dtype='>i4'), dtype='>i4')
        entry_group = wire['entry_0']
        del entry_group.attrs['awp']
        entry_group.attrs.create(
            'awp', np.asarray(6129000.0, dtype='>f8'), dtype='>f8')
        del entry_group['lang_values']
        entry_group.create_dataset(
            'lang_values', data=np.asarray([0], dtype='>i4'), dtype='>i4')
        copied = openmc.data.Reaction.from_hdf5(portable, {})
        assert copied.evaluated_product_metadata == metadata

        bad_awp = h5file.create_group('bad_awp_interpretation')
        reaction.to_hdf5(bad_awp)
        entry_group = bad_awp['evaluated_product_metadata/entry_0']
        entry_group.attrs.modify('awp_interpretation', 'arbitrary')
        with pytest.raises(ValueError):
            openmc.data.Reaction.from_hdf5(bad_awp, {})

    invalid_metadata = openmc.data.EvaluatedProductMetadata(
        102, 1, 2, 'complete', (entry,))
    invalid_reaction = openmc.data.Reaction(102)
    invalid_reaction.evaluated_product_metadata = invalid_metadata
    with h5py.File(tmp_path / 'invalid_metadata.h5', 'w') as h5file:
        with pytest.raises(ValueError):
            invalid_reaction.to_hdf5(h5file.create_group('reaction_102'))

    empty_active = openmc.data.EvaluatedProductMetadata(
        18, 1, 1, 'complete', ())
    empty_reaction = openmc.data.Reaction(18)
    empty_reaction.evaluated_product_metadata = empty_active
    with h5py.File(tmp_path / 'empty_active.h5', 'w') as h5file:
        with pytest.raises(ValueError):
            empty_reaction.to_hdf5(h5file.create_group('reaction_18'))

    mismatched_reaction = openmc.data.Reaction(51)
    mismatched_reaction.evaluated_product_metadata = metadata
    with h5py.File(tmp_path / 'mismatched_mt.h5', 'w') as h5file:
        with pytest.raises(ValueError):
            mismatched_reaction.to_hdf5(h5file.create_group('reaction_51'))

    ordered_evaluation = _metadata_section(18, 1, [
        (1, 1.0, 0, 0, ''), (1, 1.0, 0, -5, ''),
        (0, 0.0, 0, 0, ''), (0, 0.0, 0, -15, '')], jp=11)
    ordered = openmc.data.get_evaluated_product_metadata(
        ordered_evaluation, 18).metadata
    ordered_reaction = openmc.data.Reaction(18)
    ordered_reaction.evaluated_product_metadata = ordered
    with h5py.File(tmp_path / 'ordered_metadata.h5', 'w') as h5file:
        group = h5file.create_group('reaction_18')
        ordered_reaction.to_hdf5(group)
        copied = openmc.data.Reaction.from_hdf5(group, {})
        assert copied.evaluated_product_metadata == ordered
        assert tuple(entry.zap for entry in ordered.entries) == (1, 1, 0, 0)


def test_evaluated_product_metadata_law6_hdf5_wire(tmp_path):
    entry = openmc.data.EvaluatedProductMetadataEntry(
        0, 1, 1.0, 'mass_ratio', 0, 6, (), None, (), 'complete',
        'special_particle', 'neutron', 4.25, 3)
    metadata = openmc.data.EvaluatedProductMetadata(
        51, 0, 1, 'complete', (entry,))
    reaction = openmc.data.Reaction(51)
    reaction.evaluated_product_metadata = metadata
    filename = tmp_path / 'law6_metadata.h5'

    with h5py.File(filename, 'w') as h5file:
        group = h5file.create_group('law6')
        reaction.to_hdf5(group)
        entry_group = group['evaluated_product_metadata/entry_0']
        assert entry_group.attrs['apsx'] == 4.25
        assert entry_group.attrs['npsx'] == 3
        assert entry_group.attrs.get_id('apsx').dtype == np.dtype('float64')
        assert entry_group.attrs.get_id('npsx').dtype == np.dtype('int32')
        copied = openmc.data.Reaction.from_hdf5(group, {})
        assert copied.evaluated_product_metadata == metadata

    with h5py.File(filename, 'a') as h5file:
        for name, mutate in [
            ('missing_apsx', lambda attrs: attrs.__delitem__('apsx')),
            ('missing_npsx', lambda attrs: attrs.__delitem__('npsx')),
        ]:
            group = h5file.create_group(name)
            reaction.to_hdf5(group)
            attrs = group['evaluated_product_metadata/entry_0'].attrs
            mutate(attrs)
            with pytest.raises(ValueError):
                openmc.data.Reaction.from_hdf5(group, {})

        wrong_apsx_dtype = h5file.create_group('wrong_apsx_dtype')
        reaction.to_hdf5(wrong_apsx_dtype)
        attrs = wrong_apsx_dtype['evaluated_product_metadata/entry_0'].attrs
        del attrs['apsx']
        attrs.create('apsx', np.float32(4.25), dtype=np.float32)
        with pytest.raises(ValueError):
            openmc.data.Reaction.from_hdf5(wrong_apsx_dtype, {})

        wrong_apsx_rank = h5file.create_group('wrong_apsx_rank')
        reaction.to_hdf5(wrong_apsx_rank)
        attrs = wrong_apsx_rank['evaluated_product_metadata/entry_0'].attrs
        del attrs['apsx']
        attrs.create('apsx', np.asarray([4.25], dtype=np.float64),
                     dtype=np.float64)
        with pytest.raises(ValueError):
            openmc.data.Reaction.from_hdf5(wrong_apsx_rank, {})

        wrong_npsx_dtype = h5file.create_group('wrong_npsx_dtype')
        reaction.to_hdf5(wrong_npsx_dtype)
        attrs = wrong_npsx_dtype['evaluated_product_metadata/entry_0'].attrs
        del attrs['npsx']
        attrs.create('npsx', np.int64(2**31), dtype=np.int64)
        with pytest.raises(ValueError):
            openmc.data.Reaction.from_hdf5(wrong_npsx_dtype, {})

        nonfinite = h5file.create_group('nonfinite_apsx')
        reaction.to_hdf5(nonfinite)
        nonfinite['evaluated_product_metadata/entry_0'].attrs.modify(
            'apsx', np.float64(float('nan')))
        with pytest.raises(ValueError):
            openmc.data.Reaction.from_hdf5(nonfinite, {})

        nonpositive = h5file.create_group('nonpositive_apsx')
        reaction.to_hdf5(nonpositive)
        nonpositive['evaluated_product_metadata/entry_0'].attrs.modify(
            'apsx', np.float64(0.0))
        with pytest.raises(ValueError):
            openmc.data.Reaction.from_hdf5(nonpositive, {})

        small_npsx = h5file.create_group('small_npsx')
        reaction.to_hdf5(small_npsx)
        small_npsx['evaluated_product_metadata/entry_0'].attrs.modify(
            'npsx', np.int32(2))
        with pytest.raises(ValueError):
            openmc.data.Reaction.from_hdf5(small_npsx, {})

        nonlaw = h5file.create_group('nonlaw_has_law6_field')
        nonlaw_entry = openmc.data.EvaluatedProductMetadataEntry(
            0, 1, 1.0, 'mass_ratio', 0, 0, (), None, (), 'complete',
            'special_particle', 'neutron')
        nonlaw_reaction = openmc.data.Reaction(51)
        nonlaw_reaction.evaluated_product_metadata = (
            openmc.data.EvaluatedProductMetadata(
                51, 0, 1, 'complete', (nonlaw_entry,)))
        nonlaw_reaction.to_hdf5(nonlaw)
        nonlaw['evaluated_product_metadata/entry_0'].attrs.create(
            'apsx', np.float64(4.25), dtype=np.float64)
        with pytest.raises(ValueError):
            openmc.data.Reaction.from_hdf5(nonlaw, {})


def test_evaluated_product_metadata_default_is_off():
    reaction = openmc.data.Reaction(51)
    assert reaction.evaluated_product_metadata is None
    assert reaction.products == []
    parameter = inspect.signature(
        openmc.data.IncidentNeutron.from_njoy).parameters[
            'include_evaluated_product_metadata']
    assert parameter.default is False


def test_evaluated_product_metadata_from_njoy_opt_in(monkeypatch):
    class FakeLibrary:
        def __init__(self, filename):
            self.tables = [object()]

    created = []

    def fake_from_ace(cls, table):
        data = type('FakeIncidentNeutron', (), {})()
        data.reactions = {51: openmc.data.Reaction(51)}
        data.energy = {'0K': np.array([0.0])}
        created.append(data)
        return data

    monkeypatch.setattr(neutron_module, 'make_ace', lambda *args, **kwargs: None)
    monkeypatch.setattr(neutron_module, 'Library', FakeLibrary)
    monkeypatch.setattr(
        openmc.data.IncidentNeutron, 'from_ace', classmethod(fake_from_ace))

    evaluation = type('FakeEvaluation', (), {
        'gnds_name': 'Fe56', 'section': {}})()
    calls = []

    def fake_metadata(ev, mt):
        calls.append((ev, mt))
        entry = openmc.data.EvaluatedProductMetadataEntry(
            0, 1, 1.0, 'mass_ratio', 0, 0, (), None, (), 'complete',
            'special_particle', 'neutron')
        return openmc.data.EvaluatedProductMetadataResult(
            mt=mt, result_status='published_structurally_complete',
            jp=0, lct=1, section_semantic_status='complete',
            hdf5_group_version=1, entries=(entry,))

    monkeypatch.setattr(
        neutron_module, 'get_evaluated_product_metadata', fake_metadata)
    default_data = openmc.data.IncidentNeutron.from_njoy(
        'unused.endf', evaluation=evaluation, heatr=False)
    assert calls == []
    assert default_data.reactions[51].evaluated_product_metadata is None

    opted_in = openmc.data.IncidentNeutron.from_njoy(
        'unused.endf', evaluation=evaluation, heatr=False,
        include_evaluated_product_metadata=True)
    assert calls == [(evaluation, 51)]
    assert opted_in.reactions[51].evaluated_product_metadata.mt == 51
    assert opted_in.reactions[51].products == []

    monkeypatch.setattr(
        neutron_module, 'get_evaluated_product_metadata',
        lambda ev, mt: openmc.data.EvaluatedProductMetadataResult(
            mt=mt, result_status='structural_failure',
            structural_failure='unknown_law'))
    failed = openmc.data.IncidentNeutron.from_njoy(
        'unused.endf', evaluation=evaluation, heatr=False,
        include_evaluated_product_metadata=True)
    assert failed.reactions[51].evaluated_product_metadata is None
