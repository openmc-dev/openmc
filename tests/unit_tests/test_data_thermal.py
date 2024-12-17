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
    return openmc.data.ThermalScattering.from_endf(filename, nonstandard_endf=True)


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
    return openmc.data.ThermalScattering.from_endf(filename, nonstandard_endf=True)


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


@pytest.fixture
def fake_mixed_elastic():
    fake_tsl = openmc.data.ThermalScattering("c_D_in_7LiD", 1.9968, 4.9, [0.0253])
    fake_tsl.nuclides = ['H2']

    # Create elastic reaction
    bragg_edges = [0.00370672, 0.00494229, 0.00988458, 0.01359131, 0.01482688,
                   0.01976918, 0.02347589, 0.02471147, 0.02965376, 0.03336048,
                   0.03953834, 0.04324506, 0.04448063, 0.04942292, 0.05312964,
                   0.05436522, 0.05930751, 0.06301423, 0.0642498 , 0.06919209,
                   0.07289881, 0.07907667, 0.08278339, 0.08401896, 0.08896126,
                   0.09266798, 0.09390355, 0.09884584, 0.1025526 , 0.1037882 ,
                   0.1087305 , 0.1124372 , 0.1186151 , 0.1223218 , 0.1235574 ,
                   0.1284997 , 0.1322064 , 0.133442  , 0.142091  , 0.1433266 ,
                   0.1482688 , 0.1519756 , 0.1581534 , 0.1618601 , 0.1630957 ,
                   0.168038  , 0.1717447 , 0.1729803 , 0.1779226 , 0.1816293 ,
                   0.1828649 , 0.1878072 , 0.1915139 , 0.1976918 , 0.2026341 ,
                   0.2075763 , 0.2125186 , 0.2174609 , 0.2224032 , 0.2273455 ,
                   0.2421724 , 0.2471147 , 0.252057  , 0.2569993 , 0.2619415 ,
                   0.2668838 , 0.2767684 , 0.2817107 , 0.2915953 , 0.3064222 ,
                   0.3261913 , 0.366965]
    factors = [0.00375735, 0.01386287, 0.02595574, 0.02992438, 0.03549502,
               0.03855745, 0.04058831, 0.04986305, 0.05703106, 0.05855471,
               0.06078031, 0.06212291, 0.06656602, 0.06930339, 0.0697072 ,
               0.07201456, 0.07263853, 0.07313129, 0.07465531, 0.07714482,
               0.07759976, 0.077809  , 0.07790282, 0.07927957, 0.08013058,
               0.08026637, 0.08073475, 0.08112202, 0.08123039, 0.08187171,
               0.08213756, 0.08218236, 0.08236572, 0.08240729, 0.08259795,
               0.08297893, 0.08300455, 0.08314566, 0.08315611, 0.08337715,
               0.08350026, 0.08350663, 0.08352815, 0.08353776, 0.0836098 ,
               0.08367017, 0.08367361, 0.0837242 , 0.08375069, 0.08375227,
               0.08377006, 0.08381488, 0.08381644, 0.08382698, 0.08386266,
               0.08387756, 0.08388445, 0.08388974, 0.08390341, 0.08391088,
               0.08391695, 0.08392361, 0.08392684, 0.08392818, 0.08393161,
               0.08393546, 0.08393685, 0.08393801, 0.08393976, 0.08394167,
               0.08394288, 0.08394398]
    coherent_xs = openmc.data.CoherentElastic(bragg_edges, factors)
    incoherent_xs = openmc.data.Tabulated1D([0.00370672, 0.00370672], [0.00370672, 0.00370672])
    elastic_xs = {'294K': openmc.data.Sum((coherent_xs, incoherent_xs))}
    coherent_dist = openmc.data.CoherentElasticAE(coherent_xs)
    incoherent_dist = openmc.data.IncoherentElasticAEDiscrete([
        [-0.6, -0.18, 0.18, 0.6], [-0.6, -0.18, 0.18, 0.6]
    ])
    elastic_dist = {'294K': openmc.data.MixedElasticAE(coherent_dist, incoherent_dist)}
    fake_tsl.elastic = openmc.data.ThermalScatteringReaction(elastic_xs, elastic_dist)

    # Create inelastic reaction
    inelastic_xs = {'294K': openmc.data.Tabulated1D([1.0e-5, 4.9], [13.4, 3.35])}
    breakpoints = [3]
    interpolation = [2]
    energy = [1.0e-5, 4.3e-2, 4.9]
    energy_out = [
        openmc.data.Tabular([0.0002, 0.067, 0.146, 0.366], [0.25, 0.25, 0.25, 0.25]),
        openmc.data.Tabular([0.0001, 0.009, 0.137, 0.277], [0.25, 0.25, 0.25, 0.25]),
        openmc.data.Tabular([0.0579, 4.555, 4.803, 4.874], [0.25, 0.25, 0.25, 0.25]),
    ]
    for eout in energy_out:
        eout.normalize()
        eout.c = eout.cdf()
    discrete = openmc.stats.Discrete([-0.9, -0.6, -0.3, -0.1, 0.1, 0.3, 0.6, 0.9], [1/8]*8)
    discrete.c = discrete.cdf()[1:]
    mu = [[discrete]*4]*3
    inelastic_dist = {'294K': openmc.data.IncoherentInelasticAE(
        breakpoints, interpolation, energy, energy_out, mu)}
    inelastic = openmc.data.ThermalScatteringReaction(inelastic_xs, inelastic_dist)
    fake_tsl.inelastic = inelastic

    return fake_tsl


def test_mixed_elastic(fake_mixed_elastic, run_in_tmpdir):
    # Write data to HDF5 and then read back
    original = fake_mixed_elastic
    original.export_to_hdf5('c_D_in_7LiD.h5')
    copy = openmc.data.ThermalScattering.from_hdf5('c_D_in_7LiD.h5')

    # Make sure data did not change as a result of HDF5 writing/reading
    assert original == copy

    # Create modified cross_sections.xml file that includes the above data
    xs = openmc.data.DataLibrary.from_xml()
    xs.register_file('c_D_in_7LiD.h5')
    xs.export_to_xml('cross_sections_mixed.xml')

    # Create a minimal model that includes the new data and run it
    mat = openmc.Material()
    mat.add_nuclide('H2', 1.0)
    mat.add_nuclide('Li7', 1.0)
    mat.set_density('g/cm3', 1.0)
    mat.add_s_alpha_beta('c_D_in_7LiD')
    sph = openmc.Sphere(r=10.0, boundary_type="vacuum")
    cell = openmc.Cell(fill=mat, region=-sph)
    model = openmc.Model()
    model.geometry = openmc.Geometry([cell])
    model.materials = openmc.Materials([mat])
    model.materials.cross_sections = "cross_sections_mixed.xml"
    model.settings.particles = 1000
    model.settings.batches = 10
    model.settings.run_mode = 'fixed source'
    model.settings.source = openmc.IndependentSource(
        energy=openmc.stats.Discrete([3.0], [1.0])  # 3 eV source
    )
    model.run()
