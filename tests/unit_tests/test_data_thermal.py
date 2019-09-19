from collections.abc import Callable
from math import exp
import os
import random

import numpy as np
import pytest
import openmc.data

from . import needs_njoy


@pytest.fixture(scope='module')
def h2o():
    """H in H2O thermal scattering data."""
    directory = os.path.dirname(os.environ['OPENMC_CROSS_SECTIONS'])
    filename = os.path.join(directory, 'c_H_in_H2O.h5')
    return openmc.data.ThermalScattering.from_hdf5(filename)


@pytest.fixture(scope='module')
def graphite():
    """Graphite thermal scattering data."""
    directory = os.path.dirname(os.environ['OPENMC_CROSS_SECTIONS'])
    filename = os.path.join(directory, 'c_Graphite.h5')
    return openmc.data.ThermalScattering.from_hdf5(filename)


@pytest.fixture(scope='module')
def h2o_njoy():
    """H in H2O generated using NJOY."""
    endf_data = os.environ['OPENMC_ENDF_DATA']
    path_h1 = os.path.join(endf_data, 'neutrons', 'n-001_H_001.endf')
    path_h2o = os.path.join(endf_data, 'thermal_scatt', 'tsl-HinH2O.endf')
    return openmc.data.ThermalScattering.from_njoy(
        path_h1, path_h2o, temperatures=[293.6, 500.0])


@pytest.fixture(scope='module')
def hzrh():
    """H in ZrH thermal scattering data."""
    endf_data = os.environ['OPENMC_ENDF_DATA']
    filename = os.path.join(endf_data, 'thermal_scatt', 'tsl-HinZrH.endf')
    return openmc.data.ThermalScattering.from_endf(filename)


@pytest.fixture(scope='module')
def hzrh_njoy():
    """H in ZrH generated using NJOY."""
    endf_data = os.environ['OPENMC_ENDF_DATA']
    path_h1 = os.path.join(endf_data, 'neutrons', 'n-001_H_001.endf')
    path_hzrh = os.path.join(endf_data, 'thermal_scatt', 'tsl-HinZrH.endf')
    with_endf_data = openmc.data.ThermalScattering.from_njoy(
        path_h1, path_hzrh, temperatures=[296.0], iwt=0
    )
    without_endf_data = openmc.data.ThermalScattering.from_njoy(
        path_h1, path_hzrh, temperatures=[296.0], use_endf_data=False, iwt=1
    )
    return with_endf_data, without_endf_data


@pytest.fixture(scope='module')
def sio2():
    """SiO2 thermal scattering data."""
    endf_data = os.environ['OPENMC_ENDF_DATA']
    filename = os.path.join(endf_data, 'thermal_scatt', 'tsl-SiO2.endf')
    return openmc.data.ThermalScattering.from_endf(filename)


def test_h2o_attributes(h2o):
    assert h2o.name == 'c_H_in_H2O'
    assert h2o.nuclides == ['H1']
    assert h2o.temperatures == ['294K']
    assert h2o.atomic_weight_ratio == pytest.approx(0.999167)
    assert h2o.energy_max == pytest.approx(4.46)
    assert isinstance(repr(h2o), str)


def test_h2o_xs(h2o):
    assert not h2o.elastic
    for temperature, func in h2o.inelastic.xs.items():
        assert temperature.endswith('K')
        assert isinstance(func, Callable)


def test_graphite_attributes(graphite):
    assert graphite.name == 'c_Graphite'
    assert graphite.nuclides == ['C0', 'C12', 'C13']
    assert graphite.temperatures == ['296K']
    assert graphite.atomic_weight_ratio == pytest.approx(11.898)
    assert graphite.energy_max == pytest.approx(4.46)


def test_graphite_xs(graphite):
    for temperature, func in graphite.elastic.xs.items():
        assert temperature.endswith('K')
        assert isinstance(func, openmc.data.CoherentElastic)
    for temperature, func in graphite.inelastic.xs.items():
        assert temperature.endswith('K')
        assert isinstance(func, Callable)
    elastic = graphite.elastic.xs['296K']
    assert elastic([1e-3, 1.0]) == pytest.approx([0.0, 0.62586153])

@needs_njoy
def test_graphite_njoy():
    endf_data = os.environ['OPENMC_ENDF_DATA']
    path_c0 = os.path.join(endf_data, 'neutrons', 'n-006_C_000.endf')
    path_gr = os.path.join(endf_data, 'thermal_scatt', 'tsl-graphite.endf')
    graphite = openmc.data.ThermalScattering.from_njoy(
        path_c0, path_gr, temperatures=[296.0])
    assert graphite.nuclides == ['C0', 'C12', 'C13']
    assert graphite.atomic_weight_ratio == pytest.approx(11.898)
    assert graphite.energy_max == pytest.approx(2.02)
    assert graphite.temperatures == ['296K']


@needs_njoy
def test_export_to_hdf5(tmpdir, h2o_njoy, hzrh_njoy, graphite):
    filename = str(tmpdir.join('water.h5'))
    h2o_njoy.export_to_hdf5(filename)
    assert os.path.exists(filename)

    # Graphite covers export of coherent elastic data
    filename = str(tmpdir.join('graphite.h5'))
    graphite.export_to_hdf5(filename)
    assert os.path.exists(filename)

    # H in ZrH covers export of incoherent elastic data, and incoherent
    # inelastic angle-energy distributions
    filename = str(tmpdir.join('hzrh.h5'))
    hzrh_njoy[0].export_to_hdf5(filename)
    assert os.path.exists(filename)
    hzrh_njoy[1].export_to_hdf5(filename, 'w')
    assert os.path.exists(filename)


@needs_njoy
def test_continuous_dist(h2o_njoy):
    for temperature, dist in h2o_njoy.inelastic.distribution.items():
        assert temperature.endswith('K')
        assert isinstance(dist, openmc.data.IncoherentInelasticAE)


def test_h2o_endf():
    endf_data = os.environ['OPENMC_ENDF_DATA']
    filename = os.path.join(endf_data, 'thermal_scatt', 'tsl-HinH2O.endf')
    h2o = openmc.data.ThermalScattering.from_endf(filename)
    assert not h2o.elastic
    assert h2o.atomic_weight_ratio == pytest.approx(0.99917)
    assert h2o.energy_max == pytest.approx(3.99993)
    assert h2o.temperatures == ['294K', '350K', '400K', '450K', '500K', '550K',
                                '600K', '650K', '800K']


def test_hzrh_attributes(hzrh):
    assert hzrh.atomic_weight_ratio == pytest.approx(0.99917)
    assert hzrh.energy_max == pytest.approx(1.9734)
    assert hzrh.temperatures == ['296K', '400K', '500K', '600K', '700K', '800K',
                                 '1000K', '1200K']


def test_hzrh_elastic(hzrh):
    rx = hzrh.elastic
    for temperature, func in rx.xs.items():
        assert temperature.endswith('K')
        assert isinstance(func, openmc.data.IncoherentElastic)

    xs = rx.xs['296K']
    sig_b, W = xs.bound_xs, xs.debye_waller
    assert sig_b == pytest.approx(81.98006)
    assert W == pytest.approx(8.486993)
    for i in range(10):
        E = random.uniform(0.0, hzrh.energy_max)
        assert xs(E) == pytest.approx(sig_b/2 * ((1 - exp(-4*E*W))/(2*E*W)))

    for temperature, dist in rx.distribution.items():
        assert temperature.endswith('K')
        assert dist.debye_waller > 0.0


@needs_njoy
def test_hzrh_njoy(hzrh_njoy):
    endf, ace = hzrh_njoy

    # First check version using ENDF incoherent elastic data
    assert endf.atomic_weight_ratio == pytest.approx(0.999167)
    assert endf.energy_max == pytest.approx(1.855)
    assert endf.temperatures == ['296K']

    # Now check version using ACE incoherent elastic data (discretized)
    assert ace.atomic_weight_ratio == endf.atomic_weight_ratio
    assert ace.energy_max == endf.energy_max

    # Cross sections should be about the same (within 1%)
    E = np.linspace(1e-5, endf.energy_max)
    xs1 = endf.elastic.xs['296K'](E)
    xs2 = ace.elastic.xs['296K'](E)
    assert xs1 == pytest.approx(xs2, rel=0.01)

    # Check discrete incoherent elastic distribution
    d = ace.elastic.distribution['296K']
    assert np.all((-1.0 <= d.mu_out) & (d.mu_out <= 1.0))

    # Check discrete incoherent inelastic distribution
    d = endf.inelastic.distribution['296K']
    assert d.skewed
    assert np.all((-1.0 <= d.mu_out) & (d.mu_out <= 1.0))
    assert np.all((0.0 <= d.energy_out) & (d.energy_out < 3*endf.energy_max))


def test_sio2_attributes(sio2):
    assert sio2.atomic_weight_ratio == pytest.approx(27.84423)
    assert sio2.energy_max == pytest.approx(2.46675)
    assert sio2.temperatures == ['294K', '350K', '400K', '500K', '800K',
                                 '1000K', '1200K']


def test_sio2_elastic(sio2):
    rx = sio2.elastic
    for temperature, func in rx.xs.items():
        assert temperature.endswith('K')
        assert isinstance(func, openmc.data.CoherentElastic)
    xs = rx.xs['294K']
    assert len(xs) == 317
    assert xs.bragg_edges[0] == pytest.approx(0.000711634)
    assert xs.factors[0] == pytest.approx(2.6958e-14)

    # Below first bragg edge, cross section should be zero
    E = xs.bragg_edges[0] / 2.0
    assert xs(E) == 0.0

    # Between bragg edges, cross section is P/E where P is the factor
    E = (xs.bragg_edges[0] + xs.bragg_edges[1]) / 2.0
    P = xs.factors[0]
    assert xs(E) == pytest.approx(P / E)

    # Check the last Bragg edge
    E = 1.1 * xs.bragg_edges[-1]
    P = xs.factors[-1]
    assert xs(E) == pytest.approx(P / E)

    for temperature, dist in rx.distribution.items():
        assert temperature.endswith('K')
        assert dist.coherent_xs is rx.xs[temperature]


def test_get_thermal_name():
    f = openmc.data.get_thermal_name
    # Names which are recognized
    assert f('lwtr') == 'c_H_in_H2O'
    assert f('hh2o') == 'c_H_in_H2O'

    with pytest.warns(UserWarning, match='is not recognized'):
        # Names which can be guessed
        assert f('lw00') == 'c_H_in_H2O'
        assert f('graphite') == 'c_Graphite'
        assert f('D_in_D2O') == 'c_D_in_D2O'

        # Not in values, but very close
        assert f('hluci') == 'c_H_in_C5O2H8'
        assert f('ortho_d') == 'c_ortho_D'

        # Names that don't remotely match anything
        assert f('boogie_monster') == 'c_boogie_monster'
