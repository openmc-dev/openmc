#include <iostream>
#include <map>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "openmc/interpolate.h"
#include "openmc/search.h"

using namespace openmc;

TEST_CASE("Test Lagranian Interpolation")
{
  std::vector<double> xs {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  std::vector<double> ys {0.0, 1.0, 1.0, 2.0, 3.0, 3.0, 5.0};

  // ensure we get data points back at the x values
  for (int n = 1; n <= 6; n++) {
    for (int i = 0; i < xs.size(); i++) {
      double x = xs[i];
      double y = ys[i];

      size_t idx = lower_bound_index(xs.begin(), xs.end(), x);
      idx = std::min(idx, xs.size() - n - 1);
      double out = interpolate_lagrangian(xs, ys, idx, x, n);
      REQUIRE(out == y);
    }
  }

  // spot checks based on an independent implementation of Lagrangian
  // interpolation
  std::map<int, std::vector<std::pair<double, double>>> checks;
  checks[1] = {{0.5, 0.5}, {4.5, 3.0}, {2.5, 1.5}, {5.5, 4.0}};
  checks[2] = {{2.5, 1.5}, {4.5, 2.75}, {4.9999, 3.0}, {4.00001, 3.0}};
  checks[3] = {{2.5592, 1.5}, {4.5, 2.9375}, {4.9999, 3.0}, {4.00001, 3.0}};

  for (auto check_set : checks) {
    int order = check_set.first;
    auto checks = check_set.second;

    for (auto check : checks) {
      double input = check.first;
      double exp_output = check.second;

      size_t idx = lower_bound_index(xs.begin(), xs.end(), input);
      idx = std::min(idx, xs.size() - order - 1);
      double out = interpolate_lagrangian(xs, ys, idx, input, order);
      REQUIRE_THAT(out, Catch::Matchers::WithinAbs(exp_output, 1e-04));
    }
  }
}