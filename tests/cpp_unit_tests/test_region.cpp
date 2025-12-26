#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>

#include "openmc/cell.h"

TEST_CASE("Test region simplification")
{
  auto region = openmc::Region::Region("-1 2 (-3 4) | (-5 6)", 0);
  REQUIRE_THAT(region->str() == "(-1 2 -3 4) | (-5 6)");
}
