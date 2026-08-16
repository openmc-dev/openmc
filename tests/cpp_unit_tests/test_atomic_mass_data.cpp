#include <string>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "openmc/atomic_mass_data.h"

using namespace openmc::data;

TEST_CASE("Test read and storing atomic data")
{
  AtomicData my_atomic_data = AtomicData("mass_1.mas20.txt");

  std::string nucs[3] = {"Cl50", "Rh118", "Nd151"};
  double expected_mass[3] = {50.008266, 117.930341116, 150.923839363};
  double expected_binding[3] = {7651, 8322.8539, 8230.2741};
  double expected_mass_excess[3] = {7700, -64886.840, -70943.183};

  for (int i = 0; i < 3; i++) {
    double mass_out = my_atomic_data.get_atomic_mass(nucs[i]);
    REQUIRE_THAT(mass_out, Catch::Matchers::WithinAbs(expected_mass[i], 1e-04));

    double binding_out = my_atomic_data.get_atomic_binding_energy(nucs[i]);
    REQUIRE_THAT(
      binding_out, Catch::Matchers::WithinAbs(expected_binding[i], 1e-04));

    double excess_out = my_atomic_data.get_atomic_mass_excess(nucs[i]);
    REQUIRE_THAT(
      excess_out, Catch::Matchers::WithinAbs(expected_mass_excess[i], 1e-04));
  }
}