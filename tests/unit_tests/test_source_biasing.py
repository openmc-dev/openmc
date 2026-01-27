"""Tests for source biasing using C++ sampling routines via openmc.lib

This test module validates that the C++ distribution sampling implementations
correctly handle both unbiased and biased sampling when used in source
definitions. Each test:

1. Creates a minimal model with a source using a specific energy distribution
2. Uses model.sample_external_source() to generate samples via openmc.lib
3. Extracts energies from the returned particle list
4. Validates that:
   - Unbiased sampling produces the expected mean
   - Biased sampling with importance weighting produces the expected mean
   - Weights are correctly applied (non-unity for biased case)

These tests complement the Python-level tests in test_stats.py by exercising
the full C++ sampling codepath that is used during actual simulations.
"""

import numpy as np
import pytest
import openmc

from tests.unit_tests import assert_sample_mean


@pytest.fixture
def model():
    """Create a minimal model for source sampling tests."""
    sphere = openmc.Sphere(r=100.0, boundary_type='vacuum')
    cell = openmc.Cell(region=-sphere)
    geometry = openmc.Geometry([cell])
    settings = openmc.Settings(particles=100, batches=1)
    space = openmc.stats.Point()
    angle = openmc.stats.Monodirectional((1.0, 0.0, 0.0))
    settings.source = openmc.IndependentSource(space=space, angle=angle)
    return openmc.Model(geometry=geometry, settings=settings)


@pytest.mark.flaky(reruns=1)
def test_discrete(run_in_tmpdir, model):
    """Test Discrete distribution sampling via C++ routines."""
    vals = np.array([1.0, 2.0, 3.0])
    probs = np.array([0.1, 0.7, 0.2])
    exp_mean = (vals * probs).sum()

    # Create source with discrete energy distribution
    model.settings.source[0].energy = energy_dist = openmc.stats.Discrete(vals, probs)

    # Sample using C++ routines and extract energies
    n_samples = 10_000
    particles = model.sample_external_source(n_samples)
    energies = np.array([p.E for p in particles])

    # Check unbiased mean
    assert_sample_mean(energies, exp_mean)

    # Sample from biased distribution
    energy_dist.bias = np.array([0.2, 0.1, 0.7])
    particles = model.sample_external_source(n_samples)
    energies = np.array([p.E for p in particles])
    weights = np.array([p.wgt for p in particles])

    # Check biased weighted mean
    weighted_energies = energies * weights
    assert_sample_mean(weighted_energies, exp_mean)
    assert np.any(weights != 1.0)


@pytest.mark.flaky(reruns=1)
def test_uniform(run_in_tmpdir, model):
    """Test Uniform distribution sampling via C++ routines."""
    a, b = 5.0, 10.0
    exp_mean = 0.5 * (a + b)

    # Create source with uniform energy distribution
    model.settings.source[0].energy = energy_dist = openmc.stats.Uniform(a, b)

    # Sample using C++ routines and extract energies
    n_samples = 10_000
    particles = model.sample_external_source(n_samples)
    energies = np.array([p.E for p in particles])

    # Check unbiased mean
    assert_sample_mean(energies, exp_mean)

    # Sample from biased distribution
    energy_dist.bias = openmc.stats.PowerLaw(a, b, 2)
    particles = model.sample_external_source(n_samples)
    energies = np.array([p.E for p in particles])
    weights = np.array([p.wgt for p in particles])

    # Check biased weighted mean
    weighted_energies = energies * weights
    assert_sample_mean(weighted_energies, exp_mean)
    assert np.any(weights != 1.0)


@pytest.mark.flaky(reruns=1)
def test_powerlaw(run_in_tmpdir, model):
    """Test PowerLaw distribution sampling via C++ routines."""
    a, b, n = 1.0, 20.0, 2.0

    # Determine mean of distribution
    exp_mean = (n+1)*(b**(n+2) - a**(n+2))/((n+2)*(b**(n+1) - a**(n+1)))

    # Create source with powerlaw energy distribution
    model.settings.source[0].energy = energy_dist = openmc.stats.PowerLaw(a, b, n)

    # Sample using C++ routines and extract energies
    n_samples = 10_000
    particles = model.sample_external_source(n_samples)
    energies = np.array([p.E for p in particles])

    # Check unbiased mean
    assert_sample_mean(energies, exp_mean)

    # Sample from biased distribution
    energy_dist.bias = openmc.stats.Uniform(a, b)
    particles = model.sample_external_source(n_samples)
    energies = np.array([p.E for p in particles])
    weights = np.array([p.wgt for p in particles])

    # Check biased weighted mean
    weighted_energies = energies * weights
    assert_sample_mean(weighted_energies, exp_mean)
    assert np.any(weights != 1.0)


@pytest.mark.flaky(reruns=1)
def test_maxwell(run_in_tmpdir, model):
    """Test Maxwell distribution sampling via C++ routines."""
    theta = 1.2895e6
    exp_mean = 3/2 * theta

    # Create source with Maxwell energy distribution
    model.settings.source[0].energy = energy_dist = openmc.stats.Maxwell(theta)

    # Sample using C++ routines and extract energies
    n_samples = 10_000
    particles = model.sample_external_source(n_samples)
    energies = np.array([p.E for p in particles])

    # Check unbiased mean
    assert_sample_mean(energies, exp_mean)

    # Sample from biased distribution
    energy_dist.bias = openmc.stats.Maxwell(theta * 1.1)
    particles = model.sample_external_source(n_samples)
    energies = np.array([p.E for p in particles])
    weights = np.array([p.wgt for p in particles])

    # Check biased weighted mean
    weighted_energies = energies * weights
    assert_sample_mean(weighted_energies, exp_mean)
    assert np.any(weights != 1.0)


@pytest.mark.flaky(reruns=1)
def test_watt(run_in_tmpdir, model):
    """Test Watt distribution sampling via C++ routines."""
    a, b = 0.965e6, 2.29e-6
    exp_mean = 3/2 * a + a**2 * b / 4

    # Create source with Watt energy distribution
    model.settings.source[0].energy = energy_dist = openmc.stats.Watt(a, b)

    # Sample using C++ routines and extract energies
    n_samples = 10_000
    particles = model.sample_external_source(n_samples)
    energies = np.array([p.E for p in particles])

    # Check unbiased mean
    assert_sample_mean(energies, exp_mean)

    # Sample from biased distribution
    energy_dist.bias = openmc.stats.Watt(a*1.05, b)
    particles = model.sample_external_source(n_samples)
    energies = np.array([p.E for p in particles])
    weights = np.array([p.wgt for p in particles])

    # Check biased weighted mean
    weighted_energies = energies * weights
    assert_sample_mean(weighted_energies, exp_mean)
    assert np.any(weights != 1.0)


@pytest.mark.flaky(reruns=1)
def test_tabular(run_in_tmpdir, model):
    """Test Tabular distribution sampling via C++ routines."""
    # Test linear-linear sampling
    x = np.array([0.0, 5.0, 7.0, 10.0])
    p = np.array([10.0, 20.0, 5.0, 6.0])

    # Create tabular distribution and normalize to get expected mean
    model.settings.source[0].energy = energy_dist = openmc.stats.Tabular(x, p, 'linear-linear')
    energy_dist.normalize()
    exp_mean = energy_dist.mean()

    # Sample using C++ routines and extract energies
    n_samples = 10_000
    particles = model.sample_external_source(n_samples)
    energies = np.array([p.E for p in particles])

    # Check unbiased mean
    assert_sample_mean(energies, exp_mean)

    # Sample from biased distribution
    energy_dist.bias = openmc.stats.Uniform(x[0], x[-1])
    particles = model.sample_external_source(n_samples)
    energies = np.array([p.E for p in particles])
    weights = np.array([p.wgt for p in particles])

    # Check biased weighted mean
    weighted_energies = energies * weights
    assert_sample_mean(weighted_energies, exp_mean)
    assert np.any(weights != 1.0)


@pytest.mark.flaky(reruns=1)
def test_mixture(run_in_tmpdir, model):
    """Test Mixture distribution sampling via C++ routines."""
    d1 = openmc.stats.Uniform(0, 5)
    d2 = openmc.stats.Uniform(3, 7)
    p = [0.5, 0.5]

    # Create mixture energy distribution
    model.settings.source[0].energy = energy_dist = openmc.stats.Mixture(p, [d1, d2])
    exp_mean = (2.5 + 5.0) / 2

    # Sample using C++ routines and extract energies
    n_samples = 10_000
    particles = model.sample_external_source(n_samples)
    energies = np.array([p.E for p in particles])

    # Check unbiased mean
    assert_sample_mean(energies, exp_mean)

    # Sample using biased sub-distribution
    energy_dist.distribution[0].bias = openmc.stats.PowerLaw(0, 5, 2)
    particles = model.sample_external_source(n_samples)
    energies = np.array([p.E for p in particles])
    weights = np.array([p.wgt for p in particles])

    # Check biased weighted mean
    weighted_energies = energies * weights
    assert_sample_mean(weighted_energies, exp_mean)
    assert np.any(weights != 1.0)


@pytest.mark.flaky(reruns=1)
def test_normal(run_in_tmpdir, model):
    """Test Normal distribution sampling via C++ routines."""
    mean_val = 25.0
    std_dev = 2.0

    # Create source with normal energy distribution
    model.settings.source[0].energy = energy_dist = openmc.stats.Normal(mean_val, std_dev)

    # Sample using C++ routines and extract energies
    n_samples = 10_000
    particles = model.sample_external_source(n_samples)
    energies = np.array([p.E for p in particles])

    # Check unbiased mean
    assert_sample_mean(energies, mean_val)

    # Sample from biased distribution
    energy_dist.bias = openmc.stats.Normal(mean_val * 1.1, std_dev)
    particles = model.sample_external_source(n_samples)
    energies = np.array([p.E for p in particles])
    weights = np.array([p.wgt for p in particles])

    # Check biased weighted mean
    weighted_energies = energies * weights
    assert_sample_mean(weighted_energies, mean_val)
    assert np.any(weights != 1.0)
