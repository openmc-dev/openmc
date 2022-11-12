#include "openmc/distribution.h"
#include "openmc/random_lcg.h"
#include <catch2/catch_test_macros.hpp>

TEST_CASE("Test alias method sampling of a discrete distribution")
{
  int n_samples = 1e6;
  double x[5] = {-1.6, 1.1, 20.3, 4.7, 0.9};
  double p[5] = {0.2, 0.1, 0.5, 0.05, 0.15};
  double samples[n_samples];

  // Initialize distribution
  openmc::Discrete dist(x, p, 5);
  uint64_t seed = openmc::init_seed(0, 0);

  // Calculate expected distribution mean
  double mean = 0.0;
  for (size_t i = 0; i < 5; i++) {
    mean += x[i] * p[i];
  }

  // Sample distribution and calculate mean
  double dist_mean = 0.0;
  for (size_t i = 0; i < n_samples; i++) {
    samples[i] = dist.sample(&seed);
    dist_mean += samples[i];
  }
  dist_mean /= n_samples;

  // Calculate standard deviation
  double std = 0.0;
  for (size_t i = 0; i < n_samples; i++) {
    std += (samples[i] - dist_mean) * (samples[i] - dist_mean);
  }
  std /= n_samples;

  // Require sampled distribution mean is within 3 standard deviations of the
  // expected mean
  REQUIRE(std::abs(dist_mean - mean) < 3 * std);
}
