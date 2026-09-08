#include <cmath>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "openmc/array.h"
#include "openmc/math_functions.h"
#include "openmc/tensor.h"

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using openmc::array;
using openmc::combine_estimates;

namespace {

//! Covariance in which no two diagonal entries match, so that the
//! three-estimate expression is used
openmc::tensor::StaticTensor2D<double, 3, 3> distinct_cov()
{
  openmc::tensor::StaticTensor2D<double, 3, 3> cov;
  cov.fill(0.0);
  cov(0, 0) = 9.0e-8;
  cov(1, 1) = 1.6e-7;
  cov(2, 2) = 4.0e-8;
  cov(0, 1) = cov(1, 0) = 1.1e-7;
  cov(0, 2) = cov(2, 0) = 5.5e-8;
  cov(1, 2) = cov(2, 1) = 7.5e-8;
  return cov;
}

//! Covariance in which the first two estimators are identical, as in
//! multi-group mode with survival biasing
openmc::tensor::StaticTensor2D<double, 3, 3> coincident_cov()
{
  openmc::tensor::StaticTensor2D<double, 3, 3> cov;
  cov.fill(0.0);
  cov(0, 0) = cov(1, 1) = cov(0, 1) = cov(1, 0) = 4.0e-8;
  cov(2, 2) = 4.1e-8;
  cov(0, 2) = cov(2, 0) = cov(1, 2) = cov(2, 1) = 3.9e-8;
  return cov;
}

//! Combination of two estimates, transcribed from Eq. 36 and Eq. 40 of
//! Urbatsch's LA-12658-MS. Those equations are written in terms of the matrix
//! S, so the sample covariance is converted with S = (n - 1) * Sigma before
//! being substituted. This is an independent statement of the same result and
//! is what pins down the standard deviation.
array<double, 2> two_estimate_reference(double e_i, double e_j, double sigma_ii,
  double sigma_jj, double sigma_ij, int64_t n)
{
  double s_ii = (n - 1) * sigma_ii;
  double s_jj = (n - 1) * sigma_jj;
  double s_ij = (n - 1) * sigma_ij;

  double f = e_i - e_j;
  double g = s_ii + s_jj - 2.0 * s_ij;

  double mean = e_i - (s_ii - s_ij) / g * f;
  double variance =
    (s_ii * s_jj - s_ij * s_ij) * (g + n * f * f) / (n * (n - 2) * g * g);
  return {mean, std::sqrt(variance)};
}

} // namespace

TEST_CASE("Test combine_estimates with three distinct estimates")
{
  auto cov = distinct_cov();
  int64_t n = 100;

  array<double, 2> result;
  REQUIRE(combine_estimates({0.980, 0.982, 0.981}, cov, n, result));
  REQUIRE(result[1] > 0.0);

  // The weights sum to one, so estimates that all agree must combine to
  // exactly that value
  array<double, 2> agreed;
  REQUIRE(combine_estimates({0.975, 0.975, 0.975}, cov, n, agreed));
  REQUIRE_THAT(agreed[0], WithinRel(0.975, 1e-12));

  // and shifting every estimate must shift the combination by the same amount
  array<double, 2> shifted;
  REQUIRE(combine_estimates({0.990, 0.992, 0.991}, cov, n, shifted));
  REQUIRE_THAT(shifted[0] - result[0], WithinAbs(0.01, 1e-12));
}

TEST_CASE("Test combine_estimates with two coincident estimates")
{
  // The first two estimates match, so the three-estimate expression is
  // singular and the two-estimate expression must be used instead
  array<double, 3> estimates {0.950000, 0.950000, 0.950001};
  auto cov = coincident_cov();
  int64_t n = 100;

  array<double, 2> result;
  REQUIRE(combine_estimates(estimates, cov, n, result));

  auto reference = two_estimate_reference(
    estimates[0], estimates[2], cov(0, 0), cov(2, 2), cov(0, 2), n);
  REQUIRE_THAT(result[0], WithinRel(reference[0], 1e-12));
  REQUIRE_THAT(result[1], WithinRel(reference[1], 1e-12));

  // With estimates this close the combination should be no worse than the
  // best of them
  REQUIRE(result[1] <= std::sqrt(cov(2, 2) / n));
}

TEST_CASE("Test combine_estimates standard deviation scales as 1/sqrt(n)")
{
  array<double, 3> estimates {0.950000, 0.950000, 0.950001};
  auto cov = coincident_cov();

  // Holding the per-realization covariance fixed, the standard deviation of
  // the mean falls as 1/sqrt(n). This is what fails if the two-estimate
  // expression is evaluated with the sample covariance where the derivation
  // calls for S, since the neglected factor carries its own dependence on n.
  array<double, 2> low, high;
  REQUIRE(combine_estimates(estimates, cov, 100, low));
  REQUIRE(combine_estimates(estimates, cov, 400, high));
  REQUIRE_THAT(low[1] / high[1], WithinRel(2.0, 0.02));
}

TEST_CASE("Test combine_estimates rejects too few realizations")
{
  auto cov = distinct_cov();
  array<double, 2> result;

  // The three-estimate expression has an n-3 term in a denominator and the
  // two-estimate expression an n-2 term, so a combination is only defined
  // above three realizations. Callers supply their own estimate below that.
  for (int64_t n : {int64_t(0), int64_t(1), int64_t(2), int64_t(3)}) {
    REQUIRE_FALSE(combine_estimates({0.980, 0.982, 0.981}, cov, n, result));
  }
  REQUIRE(combine_estimates({0.980, 0.982, 0.981}, cov, 4, result));
}
