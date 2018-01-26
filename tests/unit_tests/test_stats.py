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
    a, b = 10, 20
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
