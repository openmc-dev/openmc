from math import pi

import numpy as np
import pytest
import openmc
import openmc.stats


def test_discrete():
    x = [0.0, 1.0, 10.0]
    p = [0.3, 0.2, 0.5]
    d = openmc.stats.Discrete(x, p)
    assert d.x == x
    assert d.p == p
    assert len(d) == len(x)
    d.to_xml_element('distribution')

    # Single point
    d2 = openmc.stats.Discrete(1e6, 1.0)
    assert d2.x == [1e6]
    assert d2.p == [1.0]
    assert len(d2) == 1


def test_uniform():
    a, b = 10.0, 20.0
    d = openmc.stats.Uniform(a, b)
    assert d.a == a
    assert d.b == b
    assert len(d) == 2

    t = d.to_tabular()
    assert t.x == [a, b]
    assert t.p == [1/(b-a), 1/(b-a)]
    assert t.interpolation == 'histogram'

    d.to_xml_element('distribution')


def test_maxwell():
    theta = 1.2895e6
    d = openmc.stats.Maxwell(theta)
    assert d.theta == theta
    assert len(d) == 1
    d.to_xml_element('distribution')


def test_watt():
    a, b = 0.965e6, 2.29e-6
    d = openmc.stats.Watt(a, b)
    assert d.a == a
    assert d.b == b
    assert len(d) == 2
    d.to_xml_element('distribution')


def test_tabular():
    x = [0.0, 5.0, 7.0]
    p = [0.1, 0.2, 0.05]
    d = openmc.stats.Tabular(x, p, 'linear-linear')
    assert d.x == x
    assert d.p == p
    assert d.interpolation == 'linear-linear'
    assert len(d) == len(x)
    d.to_xml_element('distribution')


def test_legendre():
    # Pu239 elastic scattering at 100 keV
    coeffs = [1.000e+0, 1.536e-1, 1.772e-2, 5.945e-4, 3.497e-5, 1.881e-5]
    d = openmc.stats.Legendre(coeffs)
    assert d.coefficients == pytest.approx(coeffs)
    assert len(d) == len(coeffs)

    # Integrating distribution should yield one
    mu = np.linspace(-1., 1., 1000)
    assert np.trapz(d(mu), mu) == pytest.approx(1.0, rel=1e-4)

    with pytest.raises(NotImplementedError):
        d.to_xml_element('distribution')


def test_mixture():
    d1 = openmc.stats.Uniform(0, 5)
    d2 = openmc.stats.Uniform(3, 7)
    p = [0.5, 0.5]
    mix = openmc.stats.Mixture(p, [d1, d2])
    assert mix.probability == p
    assert mix.distribution == [d1, d2]
    assert len(mix) == 4

    with pytest.raises(NotImplementedError):
        mix.to_xml_element('distribution')


def test_polar_azimuthal():
    # default polar-azimuthal should be uniform in mu and phi
    d = openmc.stats.PolarAzimuthal()
    assert isinstance(d.mu, openmc.stats.Uniform)
    assert d.mu.a == -1.
    assert d.mu.b == 1.
    assert isinstance(d.phi, openmc.stats.Uniform)
    assert d.phi.a == 0.
    assert d.phi.b == 2*pi

    mu = openmc.stats.Discrete(1., 1.)
    phi = openmc.stats.Discrete(0., 1.)
    d = openmc.stats.PolarAzimuthal(mu, phi)
    assert d.mu == mu
    assert d.phi == phi

    elem = d.to_xml_element()
    assert elem.tag == 'angle'
    assert elem.attrib['type'] == 'mu-phi'
    assert elem.find('mu') is not None
    assert elem.find('phi') is not None


def test_isotropic():
    d = openmc.stats.Isotropic()
    elem = d.to_xml_element()
    assert elem.tag == 'angle'
    assert elem.attrib['type'] == 'isotropic'


def test_monodirectional():
    d = openmc.stats.Monodirectional((1., 0., 0.))
    assert d.reference_uvw == pytest.approx((1., 0., 0.))

    elem = d.to_xml_element()
    assert elem.tag == 'angle'
    assert elem.attrib['type'] == 'monodirectional'


def test_cartesian():
    x = openmc.stats.Uniform(-10., 10.)
    y = openmc.stats.Uniform(-10., 10.)
    z = openmc.stats.Uniform(0., 20.)
    d = openmc.stats.CartesianIndependent(x, y, z)
    assert d.x == x
    assert d.y == y
    assert d.z == z

    elem = d.to_xml_element()
    assert elem.tag == 'space'
    assert elem.attrib['type'] == 'cartesian'
    assert elem.find('x') is not None
    assert elem.find('y') is not None


def test_box():
    lower_left = (-10., -10., -10.)
    upper_right = (10., 10., 10.)
    d = openmc.stats.Box(lower_left, upper_right)
    assert d.lower_left == pytest.approx(lower_left)
    assert d.upper_right == pytest.approx(upper_right)
    assert not d.only_fissionable

    elem = d.to_xml_element()
    assert elem.tag == 'space'
    assert elem.attrib['type'] == 'box'
    assert elem.find('parameters') is not None

    # only fissionable parameter
    d2 = openmc.stats.Box(lower_left, upper_right, True)
    assert d2.only_fissionable
    elem = d2.to_xml_element()
    assert elem.attrib['type'] == 'fission'


def test_point():
    p = (-4., 2., 10.)
    d = openmc.stats.Point(p)
    assert d.xyz == pytest.approx(p)

    elem = d.to_xml_element()
    assert elem.tag == 'space'
    assert elem.attrib['type'] == 'point'
    assert elem.find('parameters') is not None
