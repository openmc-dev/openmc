#include "openmc/distribution.h"
#include "openmc/random_lcg.h"
#include <catch2/catch_test_macros.hpp>
#include <pugixml.hpp>

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
  REQUIRE(std::abs(dist_mean - mean) < 4 * std);
}

TEST_CASE("Test alias sampling method for pugixml constructor")
{
  // XML doc node for Discrete contructor
  pugi::xml_document doc;
  pugi::xml_node energy = doc.append_child("energy");
  pugi::xml_node parameters = energy.append_child("parameters");
  parameters.append_child(pugi::node_pcdata)
    .set_value("17140457.745328166 1.0");

  // Initialize discrete distribution and seed
  openmc::Discrete dist(energy);
  uint64_t seed = openmc::init_seed(0, 0);

  // Assertions
  REQUIRE(dist.x().size() == 1);
  REQUIRE(dist.p().size() == 1);
  REQUIRE(dist.alias().size() == 0);
  REQUIRE(dist.x()[0] == 17140457.745328166);
  REQUIRE(dist.p()[0] == 1.0);
  REQUIRE(dist.sample(&seed) == 17140457.745328166);
}
