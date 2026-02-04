#include "openmc/distribution.h"
#include "openmc/distribution_spatial.h"
#include "openmc/position.h"
#include "openmc/random_lcg.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <pugixml.hpp>

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
    auto sample = dist.sample(&seed).first;
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
  auto sample = dist.sample(&seed).first;

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

TEST_CASE("Test construction of SpatialBox with parameters")
{
  openmc::Position ll {-1, -2, -3};
  openmc::Position ur {30, 15, 5};
  openmc::SpatialBox box(ll, ur);

  REQUIRE(box.lower_left() == openmc::Position {-1, -2, -3});
  REQUIRE(box.upper_right() == openmc::Position {30, 15, 5});
  REQUIRE_FALSE(box.only_fissionable());
}

TEST_CASE("Test Normal distribution")
{
  // Test untruncated normal distribution
  openmc::Normal normal_unbounded(0.0, 1.0);

  // Check PDF at mean (should be 1/sqrt(2*pi) ≈ 0.3989)
  REQUIRE_THAT(
    normal_unbounded.evaluate(0.0), Catch::Matchers::WithinRel(0.3989, 0.001));

  // Check that it's not truncated
  REQUIRE_FALSE(normal_unbounded.is_truncated());

  // Check accessors
  REQUIRE(normal_unbounded.mean_value() == 0.0);
  REQUIRE(normal_unbounded.std_dev() == 1.0);
  REQUIRE(normal_unbounded.lower() == -openmc::INFTY);
  REQUIRE(normal_unbounded.upper() == openmc::INFTY);
}

TEST_CASE("Test truncated Normal distribution")
{
  // Create a truncated normal: mean=0, std=1, bounds=[-1, 1]
  openmc::Normal normal_truncated(0.0, 1.0, -1.0, 1.0);

  // Check that it's truncated
  REQUIRE(normal_truncated.is_truncated());

  // Check accessors
  REQUIRE(normal_truncated.lower() == -1.0);
  REQUIRE(normal_truncated.upper() == 1.0);

  // PDF should be zero outside bounds
  REQUIRE(normal_truncated.evaluate(-2.0) == 0.0);
  REQUIRE(normal_truncated.evaluate(2.0) == 0.0);

  // PDF inside bounds should be higher than untruncated (due to
  // renormalization)
  openmc::Normal normal_unbounded(0.0, 1.0);
  REQUIRE(normal_truncated.evaluate(0.0) > normal_unbounded.evaluate(0.0));

  // The truncated PDF at mean should be approximately 0.3989 / 0.6827 ≈ 0.584
  // (0.6827 is the probability mass of N(0,1) in [-1,1])
  REQUIRE_THAT(
    normal_truncated.evaluate(0.0), Catch::Matchers::WithinRel(0.584, 0.01));
}

TEST_CASE("Test truncated Normal sampling")
{
  constexpr int n_samples = 10000;
  openmc::Normal normal_truncated(0.0, 1.0, -1.0, 1.0);
  uint64_t seed = openmc::init_seed(0, 0);

  // Sample and verify all samples are within bounds
  for (int i = 0; i < n_samples; ++i) {
    auto [x, w] = normal_truncated.sample(&seed);
    REQUIRE(x >= -1.0);
    REQUIRE(x <= 1.0);
    REQUIRE(w == 1.0); // Unbiased sampling should have weight 1
  }
}

TEST_CASE("Test one-sided truncated Normal")
{
  // Test lower-bounded only (positive half-normal)
  openmc::Normal lower_bounded(0.0, 1.0, 0.0, openmc::INFTY);
  REQUIRE(lower_bounded.is_truncated());
  REQUIRE(lower_bounded.evaluate(-1.0) == 0.0);
  REQUIRE(lower_bounded.evaluate(1.0) > 0.0);

  // PDF at 0 should be approximately 2 * 0.3989 ≈ 0.798 (half-normal)
  REQUIRE_THAT(
    lower_bounded.evaluate(0.0), Catch::Matchers::WithinRel(0.798, 0.01));

  // Test upper-bounded only
  openmc::Normal upper_bounded(0.0, 1.0, -openmc::INFTY, 0.0);
  REQUIRE(upper_bounded.is_truncated());
  REQUIRE(upper_bounded.evaluate(1.0) == 0.0);
  REQUIRE(upper_bounded.evaluate(-1.0) > 0.0);
}

TEST_CASE("Test Normal XML constructor with truncation")
{
  // XML doc node for truncated Normal
  pugi::xml_document doc;
  pugi::xml_node energy = doc.append_child("energy");
  energy.append_child("type")
    .append_child(pugi::node_pcdata)
    .set_value("normal");
  energy.append_child("parameters")
    .append_child(pugi::node_pcdata)
    .set_value("1.0e6 1.0e5 0.8e6 1.2e6");

  openmc::Normal dist(energy);
  REQUIRE(dist.mean_value() == 1.0e6);
  REQUIRE(dist.std_dev() == 1.0e5);
  REQUIRE(dist.lower() == 0.8e6);
  REQUIRE(dist.upper() == 1.2e6);
  REQUIRE(dist.is_truncated());
}
