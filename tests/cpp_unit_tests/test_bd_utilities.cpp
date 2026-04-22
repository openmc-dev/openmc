#include "openmc/random_ray/bd_utilities.h"
#include "openmc/random_ray/random_ray_simulation.h"
#include "openmc/vector.h"
#include <deque>
#include <iostream>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>

using namespace openmc;

TEST_CASE("Test rhs_backwards_difference")
{
  vector<vector<double>> ref_rhs_bd_first_order {{-68.59, -36.1},
    {-108.02, -56.0}, {-134.66666666666663, -69.33333333333331},
    {-154.66666666666666, -79.33333333333331},
    {-170.66666666666666, -87.33333333333331}, {-184.0, -94.0}};

  vector<vector<double>> ref_rhs_bd_second_order {
    {-788.5999999999999, -397.9999999999999}, {-1588.0, -797.9999999999999},
    {-2321.3333333333344, -1164.6666666666667}, {-2988.0, -1497.9999999999989},
    {-3596.88888888889, -1802.4444444444416},
    {-4156.888888888892, -2082.444444444443}};

  int vector_size = 2;
  double dt = 0.1;
  // Two functions t^2 and t^3 across x = 2, 1.9, 1.8, 1.7, 1.6, 1.5, 1.4, 1.3
  vector<std::deque<double>> test_bd_vector {
    {6.859, 5.832, 4.913, 4.096, 3.375, 2.744, 2.197},
    {3.61, 3.24, 2.89, 2.56, 2.25, 1.96, 1.69}};
  vector<double> test_rhs_bd {0.0, 0.0};
  for (int i = 0; i < 6; i++) {
    int o = i + 1;
    double d0 = rhs_backwards_difference<double>(test_bd_vector[0], o, dt);
    double d1 = rhs_backwards_difference<double>(test_bd_vector[1], o, dt);
    test_rhs_bd[0] = d0;
    test_rhs_bd[1] = d1;
    REQUIRE_THAT(
      test_rhs_bd, Catch::Matchers::Approx(ref_rhs_bd_first_order[i]));
  }
  for (int i = 0; i < 6; i++) {
    int o = i + 1;
    double rhs_d0 =
      rhs_backwards_difference<double>(test_bd_vector[0], o, dt, 2);
    double rhs_d1 =
      rhs_backwards_difference<double>(test_bd_vector[1], o, dt, 2);
    test_rhs_bd[0] = rhs_d0;
    test_rhs_bd[1] = rhs_d1;
    REQUIRE_THAT(
      test_rhs_bd, Catch::Matchers::Approx(ref_rhs_bd_second_order[i]));
  }
}
