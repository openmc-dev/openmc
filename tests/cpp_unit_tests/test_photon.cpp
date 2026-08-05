#include <cmath>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "openmc/constants.h"
#include "openmc/photon.h"

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

TEST_CASE("Compton profile exponential tail")
{
  constexpr double pz_last = 5.0;
  constexpr double profile_last = 0.1;
  constexpr double slope = -0.2;

  REQUIRE(openmc::detail::compton_profile_tail_integral(
            pz_last, pz_last, profile_last, slope) == 0.0);

  double integral = openmc::detail::compton_profile_tail_integral(
    7.0, pz_last, profile_last, slope);
  double expected = profile_last * std::expm1(2.0 * slope) / slope;
  REQUIRE_THAT(integral, WithinRel(expected, 1.0e-14));

  double pz = openmc::detail::invert_compton_profile_tail(
    integral, pz_last, profile_last, slope);
  REQUIRE_THAT(pz, WithinAbs(7.0, 1.0e-13));

  double total_tail = -profile_last / slope;
  REQUIRE_THAT(openmc::detail::compton_profile_tail_integral(
                 200.0, pz_last, profile_last, slope),
    WithinRel(total_tail, 1.0e-14));
}

TEST_CASE("Compton energy root follows signed electron momentum")
{
  constexpr double alpha = 1.0;
  constexpr double mu = 0.0;
  double free_electron_ratio = 1.0 / (1.0 + alpha * (1.0 - mu));

  REQUIRE_THAT(openmc::detail::compton_energy_ratio(alpha, mu, 0.0),
    WithinAbs(free_electron_ratio, 1.0e-15));

  double negative_pz_ratio =
    openmc::detail::compton_energy_ratio(alpha, mu, -10.0);
  double positive_pz_ratio =
    openmc::detail::compton_energy_ratio(alpha, mu, 10.0);

  REQUIRE(negative_pz_ratio > 0.0);
  REQUIRE(negative_pz_ratio < free_electron_ratio);
  REQUIRE(positive_pz_ratio > free_electron_ratio);
}
