from math import pi

import numpy as np
import pytest
import openmc
import openmc.stats
from scipy.integrate import trapezoid


def assert_sample_mean(samples, expected_mean):
    # Calculate sample standard deviation
    std_dev = samples.std() / np.sqrt(samples.size - 1)

    # Means should agree within 4 sigma 99.993% of the time. Note that this is
    # expected to fail about 1 out of 16,000 times
    assert np.abs(expected_mean - samples.mean()) < 4*std_dev


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

    # sample discrete distribution and check that the mean of the samples is
    # within 4 std. dev. of the expected mean
    n_samples = 1_000_000
    samples = d3.sample(n_samples)
    assert_sample_mean(samples, exp_mean)


def test_delta_function():
    d = openmc.stats.delta_function(14.1e6)
    assert isinstance(d, openmc.stats.Discrete)
    np.testing.assert_array_equal(d.x, [14.1e6])
    np.testing.assert_array_equal(d.p, [1.0])


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
    assert merged.integral() == pytest.approx(1.0)

    # Probabilities add up but are not normalized
    d1 = openmc.stats.Discrete([3.0], [1.0])
    triple = openmc.stats.Discrete.merge([d1, d1, d1], [1.0, 2.0, 3.0])
    assert triple.x == pytest.approx([3.0])
    assert triple.p == pytest.approx([6.0])
    assert triple.integral() == pytest.approx(6.0)


def test_clip_discrete():
    # Create discrete distribution with two points that are not important, one
    # because the x value is very small, and one because the p value is very
    # small
    d = openmc.stats.Discrete([1e-8, 1.0, 2.0, 1000.0], [3.0, 2.0, 5.0, 1e-12])

    # Clipping the distribution should result in two points
    d_clip = d.clip(1e-6)
    assert d_clip.x.size == 2
    assert d_clip.p.size == 2

    # Make sure inplace returns same object
    d_same = d.clip(1e-6, inplace=True)
    assert d_same is d

    with pytest.raises(ValueError):
        d.clip(-1.)

    with pytest.raises(ValueError):
        d.clip(5)


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

    # Sample distribution and check that the mean of the samples is within 4
    # std. dev. of the expected mean
    exp_mean = 0.5 * (a + b)
    n_samples = 1_000_000
    samples = d.sample(n_samples)
    assert_sample_mean(samples, exp_mean)


def test_powerlaw():
    a, b, n = 10.0, 100.0, 2.0
    d = openmc.stats.PowerLaw(a, b, n)
    elem = d.to_xml_element('distribution')

    d = openmc.stats.PowerLaw.from_xml_element(elem)
    assert d.a == a
    assert d.b == b
    assert d.n == n
    assert len(d) == 3

    # Determine mean of distribution
    exp_mean = (n+1)*(b**(n+2) - a**(n+2))/((n+2)*(b**(n+1) - a**(n+1)))

    # sample power law distribution and check that the mean of the samples is
    # within 4 std. dev. of the expected mean
    n_samples = 1_000_000
    samples = d.sample(n_samples)
    assert_sample_mean(samples, exp_mean)


def test_maxwell():
    theta = 1.2895e6
    d = openmc.stats.Maxwell(theta)
    elem = d.to_xml_element('distribution')

    d = openmc.stats.Maxwell.from_xml_element(elem)
    assert d.theta == theta
    assert len(d) == 1

    exp_mean = 3/2 * theta

    # sample maxwell distribution and check that the mean of the samples is
    # within 4 std. dev. of the expected mean
    n_samples = 1_000_000
    samples = d.sample(n_samples)
    assert_sample_mean(samples, exp_mean)

    # A second sample starting from a different seed
    samples_2 = d.sample(n_samples)
    assert_sample_mean(samples_2, exp_mean)
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

    # sample Watt distribution and check that the mean of the samples is within
    # 4 std. dev. of the expected mean
    n_samples = 1_000_000
    samples = d.sample(n_samples)
    assert_sample_mean(samples, exp_mean)


def test_tabular():
    # test linear-linear sampling
    x = np.array([0.0, 5.0, 7.0, 10.0])
    p = np.array([10.0, 20.0, 5.0, 6.0])
    d = openmc.stats.Tabular(x, p, 'linear-linear')
    n_samples = 100_000
    samples = d.sample(n_samples)
    assert_sample_mean(samples, d.mean())

    # test linear-linear normalization
    d.normalize()
    assert d.integral() == pytest.approx(1.0)

    # test histogram sampling
    d = openmc.stats.Tabular(x, p, interpolation='histogram')
    samples = d.sample(n_samples)
    assert_sample_mean(samples, d.mean())

    d.normalize()
    assert d.integral() == pytest.approx(1.0)

    # ensure that passing a set of probabilities shorter than x works
    # for histogram interpolation
    d = openmc.stats.Tabular(x, p[:-1], interpolation='histogram')
    d.cdf()
    d.mean()
    assert_sample_mean(d.sample(n_samples), d.mean())

    # passing a shorter probability set should raise an error for linear-linear
    with pytest.raises(ValueError):
        d = openmc.stats.Tabular(x, p[:-1], interpolation='linear-linear')
        d.cdf()

    # Use probabilities of correct length for linear-linear interpolation and
    # call the CDF method
    d = openmc.stats.Tabular(x, p, interpolation='linear-linear')
    d.cdf()


def test_tabular_from_xml():
    x = np.array([0.0, 5.0, 7.0, 10.0])
    p = np.array([10.0, 20.0, 5.0, 6.0])
    d = openmc.stats.Tabular(x, p, 'linear-linear')
    elem = d.to_xml_element('distribution')

    d = openmc.stats.Tabular.from_xml_element(elem)
    assert all(d.x == x)
    assert all(d.p == p)
    assert d.interpolation == 'linear-linear'
    assert len(d) == len(x)

    # Make sure XML roundtrip works with len(x) == len(p) + 1
    x = np.array([0.0, 5.0, 7.0, 10.0])
    p = np.array([10.0, 20.0, 5.0])
    d = openmc.stats.Tabular(x, p, 'histogram')
    elem = d.to_xml_element('distribution')
    d = openmc.stats.Tabular.from_xml_element(elem)
    assert all(d.x == x)
    assert all(d.p == p)


def test_legendre():
    # Pu239 elastic scattering at 100 keV
    coeffs = [1.000e+0, 1.536e-1, 1.772e-2, 5.945e-4, 3.497e-5, 1.881e-5]
    d = openmc.stats.Legendre(coeffs)
    assert d.coefficients == pytest.approx(coeffs)
    assert len(d) == len(coeffs)

    # Integrating distribution should yield one
    mu = np.linspace(-1., 1., 1000)
    assert trapezoid(d(mu), mu) == pytest.approx(1.0, rel=1e-4)

    with pytest.raises(NotImplementedError):
        d.to_xml_element('distribution')


def test_mixture():
    d1 = openmc.stats.Uniform(0, 5)
    d2 = openmc.stats.Uniform(3, 7)
    p = [0.5, 0.5]
    mix = openmc.stats.Mixture(p, [d1, d2])
    np.testing.assert_allclose(mix.probability, p)
    assert mix.distribution == [d1, d2]
    assert len(mix) == 4

    # Sample and make sure sample mean is close to expected mean
    n_samples = 1_000_000
    samples = mix.sample(n_samples)
    assert_sample_mean(samples, (2.5 + 5.0)/2)

    elem = mix.to_xml_element('distribution')

    d = openmc.stats.Mixture.from_xml_element(elem)
    np.testing.assert_allclose(d.probability, p)
    assert d.distribution == [d1, d2]
    assert len(d) == 4


def test_mixture_clip():
    # Create mixture distribution containing a discrete distribution with two
    # points that are not important, one because the x value is very small, and
    # one because the p value is very small
    d1 = openmc.stats.Discrete([1e-8, 1.0, 2.0, 1000.0], [3.0, 2.0, 5.0, 1e-12])
    d2 = openmc.stats.Uniform(0, 5)
    mix = openmc.stats.Mixture([0.5, 0.5], [d1, d2])

    # Clipping should reduce the contained discrete distribution to 2 points
    mix_clip = mix.clip(1e-6)
    assert mix_clip.distribution[0].x.size == 2
    assert mix_clip.distribution[0].p.size == 2

    # Make sure inplace returns same object
    mix_same = mix.clip(1e-6, inplace=True)
    assert mix_same is mix

    # Make sure clip removes low probability distributions
    d_small = openmc.stats.Uniform(0., 1.)
    d_large = openmc.stats.Uniform(2., 5.)
    mix = openmc.stats.Mixture([1e-10, 1.0], [d_small, d_large])
    mix_clip = mix.clip(1e-3)
    assert mix_clip.distribution == [d_large]

    # Make sure warning is raised if tolerance is exceeded
    d1 = openmc.stats.Discrete([1.0, 1.001], [1.0, 0.7e-6])
    d2 = openmc.stats.Tabular([0.0, 1.0], [0.7e-6], interpolation='histogram')
    mix = openmc.stats.Mixture([1.0, 1.0], [d1, d2])
    with pytest.warns(UserWarning):
        mix_clip = mix.clip(1e-6)


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
    n_samples = 100_000
    samples = d.sample(n_samples)
    assert_sample_mean(samples, mean)


def test_combine_normal():
    norm1=openmc.stats.Normal(mean_value=10, std_dev=1)
    norm2=openmc.stats.Normal(mean_value=1, std_dev=1)
    combined=openmc.stats.combine_distributions([norm1,norm2],probs=[0.1,0.9])
    assert isinstance(combined, openmc.stats.Normal)
    assert combined.mean_value == pytest.approx(1.9)
    assert combined.std_dev**2 == pytest.approx(0.82)


def test_muir():
    mean = 10.0
    mass = 5.0
    temp = 20000.
    d = openmc.stats.muir(mean, mass, temp)
    assert isinstance(d, openmc.stats.Normal)

    elem = d.to_xml_element('energy')
    assert elem.attrib['type'] == 'normal'

    d = openmc.stats.Univariate.from_xml_element(elem)
    assert isinstance(d, openmc.stats.Normal)

    # sample muir distribution
    n_samples = 100_000
    samples = d.sample(n_samples)
    assert_sample_mean(samples, mean)


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

    # Combine 1 discrete and 2 tabular -- the tabular distributions should
    # combine to produce a uniform distribution with mean 0.5. The combined
    # distribution should have a mean of 0.25.
    t1 = openmc.stats.Tabular([0., 1.], [2.0, 0.0])
    t2 = openmc.stats.Tabular([0., 1.], [0.0, 2.0])
    d1 = openmc.stats.Discrete([0.0], [1.0])
    combined = openmc.stats.combine_distributions([t1, t2, d1], [0.25, 0.25, 0.5])
    assert combined.integral() == pytest.approx(1.0)

    # Sample the combined distribution and make sure the sample mean is within
    # uncertainty of the expected value
    samples = combined.sample(10_000)
    assert_sample_mean(samples, 0.25)
