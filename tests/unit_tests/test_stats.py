from math import pi

import numpy as np
import pytest
import openmc
import openmc.stats


def test_discrete():
    x = [0.0, 1.0, 10.0]
    p = [0.3, 0.2, 0.5]
    d = openmc.stats.Discrete(x, p)
    elem = d.to_xml_element('distribution')

    d = openmc.stats.Discrete.from_xml_element(elem)
    np.testing.assert_array_equal(d.x, x)
    np.testing.assert_array_equal(d.p, p)
    assert len(d) == len(x)

    d = openmc.stats.Univariate.from_xml_element(elem)
    assert isinstance(d, openmc.stats.Discrete)

    # Single point
    d2 = openmc.stats.Discrete(1e6, 1.0)
    assert d2.x == [1e6]
    assert d2.p == [1.0]
    assert len(d2) == 1

    vals = np.array([1.0, 2.0, 3.0])
    probs = np.array([0.1, 0.7, 0.2])

    exp_mean = (vals * probs).sum()

    d3 = openmc.stats.Discrete(vals, probs)

    # sample discrete distribution
    n_samples = 1_000_000
    samples = d3.sample(n_samples, seed=100)
    # check that the mean of the samples is within 3 std. dev.
    # of the expected mean
    std_dev = samples.std() / np.sqrt(n_samples)
    assert np.abs(exp_mean - samples.mean()) < 3*std_dev


def test_merge_discrete():
    x1 = [0.0, 1.0, 10.0]
    p1 = [0.3, 0.2, 0.5]
    d1 = openmc.stats.Discrete(x1, p1)

    x2 = [0.5, 1.0, 5.0]
    p2 = [0.4, 0.5, 0.1]
    d2 = openmc.stats.Discrete(x2, p2)

    # Merged distribution should have x values sorted and probabilities
    # appropriately combined. Duplicate x values should appear once.
    merged = openmc.stats.Discrete.merge([d1, d2], [0.6, 0.4])
    assert merged.x == pytest.approx([0.0, 0.5, 1.0, 5.0, 10.0])
    assert merged.p == pytest.approx(
        [0.6*0.3, 0.4*0.4, 0.6*0.2 + 0.4*0.5, 0.4*0.1, 0.6*0.5])

    # Probabilities add up but are not normalized
    d1 = openmc.stats.Discrete([3.0], [1.0])
    triple = openmc.stats.Discrete.merge([d1, d1, d1], [1.0, 2.0, 3.0])
    assert triple.x == pytest.approx([3.0])
    assert triple.p == pytest.approx([6.0])


def test_uniform():
    a, b = 10.0, 20.0
    d = openmc.stats.Uniform(a, b)
    elem = d.to_xml_element('distribution')

    d = openmc.stats.Uniform.from_xml_element(elem)
    assert d.a == a
    assert d.b == b
    assert len(d) == 2

    t = d.to_tabular()
    np.testing.assert_array_equal(t.x, [a, b])
    np.testing.assert_array_equal(t.p, [1/(b-a), 1/(b-a)])
    assert t.interpolation == 'histogram'

    exp_mean = 0.5 * (a + b)
    n_samples = 1_000_000
    samples = d.sample(n_samples, seed=100)
    # check that the mean of the samples is within 3 std. dev.
    # of the expected mean
    std_dev = samples.std() / np.sqrt(n_samples)
    assert np.abs(exp_mean - samples.mean()) < 3*std_dev


def test_powerlaw():
    a, b, n = 10.0, 100.0, 2.0
    d = openmc.stats.PowerLaw(a, b, n)
    elem = d.to_xml_element('distribution')

    d = openmc.stats.PowerLaw.from_xml_element(elem)
    assert d.a == a
    assert d.b == b
    assert d.n == n
    assert len(d) == 3

    exp_mean = 100.0 * (n+1) / (n+2)

    # sample power law distribution
    n_samples = 1_000_000
    samples = d.sample(n_samples, seed=100)
    # check that the mean of the samples is within 3 std. dev.
    # of the expected mean
    std_dev = samples.std() / np.sqrt(n_samples)
    assert np.abs(exp_mean - samples.mean()) < 3*std_dev


def test_maxwell():
    theta = 1.2895e6
    d = openmc.stats.Maxwell(theta)
    elem = d.to_xml_element('distribution')

    d = openmc.stats.Maxwell.from_xml_element(elem)
    assert d.theta == theta
    assert len(d) == 1

    exp_mean = 3/2 * theta

    # sample maxwell distribution
    n_samples = 1_000_000
    samples = d.sample(n_samples, seed=100)
    # check that the mean of the samples is within 3 std. dev.
    # of the expected mean
    std_dev = samples.std() / np.sqrt(n_samples)
    assert np.abs(exp_mean - samples.mean()) < 3*std_dev

    # A second sample with a different seed
    samples_2 = d.sample(n_samples, seed=200)
    # check that the mean of the samples is within 3 std. dev.
    # of the expected mean
    std_dev = samples_2.std() / np.sqrt(n_samples)
    assert np.abs(exp_mean - samples_2.mean()) < 3*std_dev

    assert samples_2.mean() != samples.mean()


def test_watt():
    a, b = 0.965e6, 2.29e-6
    d = openmc.stats.Watt(a, b)
    elem = d.to_xml_element('distribution')

    d = openmc.stats.Watt.from_xml_element(elem)
    assert d.a == a
    assert d.b == b
    assert len(d) == 2

    # mean value form adapted from
    # "Prompt-fission-neutron average energy for 238U(n, f ) from
    # threshold to 200 MeV" Ethvignot et. al.
    # https://doi.org/10.1016/j.physletb.2003.09.048
    exp_mean = 3/2 * a + a**2 * b / 4

    # sample Watt distribution
    n_samples = 1_000_000
    samples = d.sample(n_samples, seed=100)
    # check that the mean of the samples is within 3 std. dev.
    # of the expected mean
    std_dev = samples.std() / np.sqrt(n_samples)
    assert np.abs(exp_mean - samples.mean()) < 3*std_dev


def test_tabular():
    x = np.array([0.0, 5.0, 7.0])
    p = np.array([0.1, 0.2, 0.05])
    d = openmc.stats.Tabular(x, p, 'linear-linear')
    elem = d.to_xml_element('distribution')

    d = openmc.stats.Tabular.from_xml_element(elem)
    assert all(d.x == x)
    assert all(d.p == p)
    assert d.interpolation == 'linear-linear'
    assert len(d) == len(x)

    # test linear-linear sampling
    d = openmc.stats.Tabular(x, p)

    n_samples = 100_000
    samples = d.sample(n_samples, seed=100)
    diff = np.abs(samples - d.mean())
    # within_1_sigma = np.count_nonzero(diff < samples.std())
    # assert within_1_sigma / n_samples >= 0.68
    within_2_sigma = np.count_nonzero(diff < 2*samples.std())
    assert within_2_sigma / n_samples >= 0.95
    within_3_sigma = np.count_nonzero(diff < 3*samples.std())
    assert within_3_sigma / n_samples >= 0.99

    # test histogram sampling
    d = openmc.stats.Tabular(x, p, interpolation='histogram')
    d.normalize()

    samples = d.sample(n_samples, seed=100)
    diff = np.abs(samples - d.mean())
    # within_1_sigma = np.count_nonzero(diff < samples.std())
    # assert within_1_sigma / n_samples >= 0.68
    within_2_sigma = np.count_nonzero(diff < 2*samples.std())
    assert within_2_sigma / n_samples >= 0.95
    within_3_sigma = np.count_nonzero(diff < 3*samples.std())
    assert within_3_sigma / n_samples >= 0.99


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

    elem = mix.to_xml_element('distribution')

    d = openmc.stats.Mixture.from_xml_element(elem)
    assert d.probability == p
    assert d.distribution == [d1, d2]
    assert len(d) == 4


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

    d = openmc.stats.PolarAzimuthal.from_xml_element(elem)
    assert d.mu.x == [1.]
    assert d.mu.p == [1.]
    assert d.phi.x == [0.]
    assert d.phi.p == [1.]

    d = openmc.stats.UnitSphere.from_xml_element(elem)
    assert isinstance(d, openmc.stats.PolarAzimuthal)


def test_isotropic():
    d = openmc.stats.Isotropic()
    elem = d.to_xml_element()
    assert elem.tag == 'angle'
    assert elem.attrib['type'] == 'isotropic'

    d = openmc.stats.Isotropic.from_xml_element(elem)
    assert isinstance(d, openmc.stats.Isotropic)


def test_monodirectional():
    d = openmc.stats.Monodirectional((1., 0., 0.))
    elem = d.to_xml_element()
    assert elem.tag == 'angle'
    assert elem.attrib['type'] == 'monodirectional'

    d = openmc.stats.Monodirectional.from_xml_element(elem)
    assert d.reference_uvw == pytest.approx((1., 0., 0.))


def test_cartesian():
    x = openmc.stats.Uniform(-10., 10.)
    y = openmc.stats.Uniform(-10., 10.)
    z = openmc.stats.Uniform(0., 20.)
    d = openmc.stats.CartesianIndependent(x, y, z)

    elem = d.to_xml_element()
    assert elem.tag == 'space'
    assert elem.attrib['type'] == 'cartesian'
    assert elem.find('x') is not None
    assert elem.find('y') is not None

    d = openmc.stats.CartesianIndependent.from_xml_element(elem)
    assert d.x == x
    assert d.y == y
    assert d.z == z

    d = openmc.stats.Spatial.from_xml_element(elem)
    assert isinstance(d, openmc.stats.CartesianIndependent)


def test_box():
    lower_left = (-10., -10., -10.)
    upper_right = (10., 10., 10.)
    d = openmc.stats.Box(lower_left, upper_right)

    elem = d.to_xml_element()
    assert elem.tag == 'space'
    assert elem.attrib['type'] == 'box'
    assert elem.find('parameters') is not None

    d = openmc.stats.Box.from_xml_element(elem)
    assert d.lower_left == pytest.approx(lower_left)
    assert d.upper_right == pytest.approx(upper_right)
    assert not d.only_fissionable

    # only fissionable parameter
    d2 = openmc.stats.Box(lower_left, upper_right, True)
    assert d2.only_fissionable
    elem = d2.to_xml_element()
    assert elem.attrib['type'] == 'fission'
    d = openmc.stats.Spatial.from_xml_element(elem)
    assert isinstance(d, openmc.stats.Box)


def test_point():
    p = (-4., 2., 10.)
    d = openmc.stats.Point(p)

    elem = d.to_xml_element()
    assert elem.tag == 'space'
    assert elem.attrib['type'] == 'point'
    assert elem.find('parameters') is not None

    d = openmc.stats.Point.from_xml_element(elem)
    assert d.xyz == pytest.approx(p)


def test_normal():
    mean = 10.0
    std_dev = 2.0
    d = openmc.stats.Normal(mean,std_dev)

    elem = d.to_xml_element('distribution')
    assert elem.attrib['type'] == 'normal'

    d = openmc.stats.Normal.from_xml_element(elem)
    assert d.mean_value == pytest.approx(mean)
    assert d.std_dev == pytest.approx(std_dev)
    assert len(d) == 2

    # sample normal distribution
    n_samples = 10000
    samples = d.sample(n_samples, seed=100)
    samples = np.abs(samples - mean)
    within_1_sigma = np.count_nonzero(samples < std_dev)
    assert within_1_sigma / n_samples >= 0.68
    within_2_sigma = np.count_nonzero(samples < 2*std_dev)
    assert within_2_sigma / n_samples >= 0.95
    within_3_sigma = np.count_nonzero(samples < 3*std_dev)
    assert within_3_sigma / n_samples >= 0.99


def test_muir():
    mean = 10.0
    mass = 5.0
    temp = 20000.
    d = openmc.stats.Muir(mean,mass,temp)

    elem = d.to_xml_element('energy')
    assert elem.attrib['type'] == 'muir'

    d = openmc.stats.Muir.from_xml_element(elem)
    assert d.e0 == pytest.approx(mean)
    assert d.m_rat == pytest.approx(mass)
    assert d.kt == pytest.approx(temp)
    assert len(d) == 3

    # sample muir distribution
    n_samples = 10000
    samples = d.sample(n_samples, seed=100)
    samples = np.abs(samples - mean)
    within_1_sigma = np.count_nonzero(samples < d.std_dev)
    assert within_1_sigma / n_samples >= 0.68
    within_2_sigma = np.count_nonzero(samples < 2*d.std_dev)
    assert within_2_sigma / n_samples >= 0.95
    within_3_sigma = np.count_nonzero(samples < 3*d.std_dev)
    assert within_3_sigma / n_samples >= 0.99


def test_combine_distributions():
    # Combine two discrete (same data as in test_merge_discrete)
    x1 = [0.0, 1.0, 10.0]
    p1 = [0.3, 0.2, 0.5]
    d1 = openmc.stats.Discrete(x1, p1)
    x2 = [0.5, 1.0, 5.0]
    p2 = [0.4, 0.5, 0.1]
    d2 = openmc.stats.Discrete(x2, p2)

    # Merged distribution should have x values sorted and probabilities
    # appropriately combined. Duplicate x values should appear once.
    merged = openmc.stats.combine_distributions([d1, d2], [0.6, 0.4])
    assert isinstance(merged, openmc.stats.Discrete)
    assert merged.x == pytest.approx([0.0, 0.5, 1.0, 5.0, 10.0])
    assert merged.p == pytest.approx(
        [0.6*0.3, 0.4*0.4, 0.6*0.2 + 0.4*0.5, 0.4*0.1, 0.6*0.5])

    # Probabilities add up but are not normalized
    d1 = openmc.stats.Discrete([3.0], [1.0])
    triple = openmc.stats.combine_distributions([d1, d1, d1], [1.0, 2.0, 3.0])
    assert triple.x == pytest.approx([3.0])
    assert triple.p == pytest.approx([6.0])

    # Combine discrete and tabular
    t1 = openmc.stats.Tabular(x2, p2)
    mixed = openmc.stats.combine_distributions([d1, t1], [0.5, 0.5])
    assert isinstance(mixed, openmc.stats.Mixture)
    assert len(mixed.distribution) == 2
    assert len(mixed.probability) == 2
