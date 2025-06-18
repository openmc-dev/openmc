#include "openmc/distribution.h"
#include "openmc/random_lcg.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <pugixml.hpp>

TEST_CASE("Test integral over normalized distributions is 1.0")
{
  uint64_t seed = openmc::init_seed(0, 0);

  // Initialize discrete distribution
  double x[5] = {-1.6, 1.1, 20.3, 4.7, 0.9};
  double p[5] = {0.2, 0.1, 0.65, 0.02, 0.03};
  openmc::Discrete discrete(x, p, 5);

  REQUIRE_THAT(
    discrete.integral(-2.0, 21.0), Catch::Matchers::WithinAbs(1.0, 1e-12));

  openmc::Uniform uniform(-1.0, 1.0);
  REQUIRE_THAT(
    uniform.integral(-1.0, 1.0), Catch::Matchers::WithinAbs(1.0, 1e-12));

  double n[4] = {0.0, -1.0, -2.0, 3.0};
  for (auto i : n) {
    openmc::PowerLaw powerlaw(1.0, 5.0, i);
    REQUIRE_THAT(
      powerlaw.integral(1.0, 5.0), Catch::Matchers::WithinAbs(1.0, 1e-12));
  }

  openmc::Maxwell maxwell(1.0);
  REQUIRE_THAT(maxwell.integral(0.0, openmc::INFTY),
    Catch::Matchers::WithinAbs(1.0, 1e-12));

  openmc::Watt watt(0.5, 3.0);
  REQUIRE_THAT(
    watt.integral(0.0, openmc::INFTY), Catch::Matchers::WithinAbs(1.0, 1e-12));

  openmc::Normal normal(0.5, 3.0);
  REQUIRE_THAT(normal.integral(-openmc::INFTY, openmc::INFTY),
    Catch::Matchers::WithinAbs(1.0, 1e-12));
}

TEST_CASE("Test alias method sampling of a discrete distribution")
{
  constexpr int n_samples = 1000000;
  double x[5] = {-1.6, 1.1, 20.3, 4.7, 0.9};
  double p[5] = {0.2, 0.1, 0.65, 0.02, 0.03};

  // Initialize distribution
  openmc::Discrete dist(x, p, 5);
  uint64_t seed = openmc::init_seed(0, 0);

  // Calculate expected distribution mean
  double mean = 0.0;
  for (size_t i = 0; i < 5; i++) {
    mean += x[i] * p[i];
  }

  // Sample distribution and calculate mean, standard deviation, and number of
  // x[0] sampled
  double dist_mean = 0.0;
  double std = 0.0;
  int counter = 0;

  for (size_t i = 0; i < n_samples; i++) {
    auto sample = dist.sample(&seed);
    std += sample * sample / n_samples;
    dist_mean += sample;

    if (sample == x[0])
      counter++;
  }

  dist_mean /= n_samples;
  std -= dist_mean * dist_mean;

  // Require sampled distribution mean is within 4 standard deviations of the
  // expected mean
  REQUIRE(std::abs(dist_mean - mean) < 4 * std);

  // Require counter of number of x[0] is within the 95% confidence interval
  // assuming a Poisson distribution of 200,000
  REQUIRE(std::abs((double)counter / n_samples - p[0]) <
          1.96 * std::sqrt(p[0] / n_samples));
}

TEST_CASE("Test alias sampling method for pugixml constructor")
{
  // XML doc node for Discrete contructor
  pugi::xml_document doc;
  pugi::xml_node energy = doc.append_child("energy");
  pugi::xml_node parameters = energy.append_child("parameters");
  parameters.append_child(pugi::node_pcdata)
    .set_value("800 500000 30000 0.1 0.6 0.3");

  // Initialize discrete distribution and seed
  openmc::Discrete dist(energy);
  uint64_t seed = openmc::init_seed(0, 0);
  auto sample = dist.sample(&seed);

  // Assertions
  REQUIRE(dist.x().size() == 3);
  REQUIRE(dist.prob().size() == 3);
  REQUIRE(dist.alias().size() == 3);

  openmc::vector<double> correct_x = {800, 500000, 30000};
  openmc::vector<double> correct_prob = {0.3, 1.0, 0.9};
  openmc::vector<size_t> correct_alias = {1, 0, 1};

  for (size_t i = 0; i < 3; i++) {
    REQUIRE(dist.x()[i] == correct_x[i]);
    REQUIRE_THAT(
      dist.prob()[i], Catch::Matchers::WithinAbs(correct_prob[i], 1e-12));
    REQUIRE(dist.alias()[i] == correct_alias[i]);
  }
}
