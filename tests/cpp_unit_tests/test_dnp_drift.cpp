#include <catch2/catch_test_macros.hpp>

#include "openmc/dnp_drift.h"
#include "openmc/position.h"

using namespace openmc;

TEST_CASE("Test linear position adjustment")
{
  // Typical use case
  Position a = Position(-1.0, -1.0, -1.0);
  Position b = Position(1.0, 1.0, 1.0);
  _adjust_position(a, b, 3.0, 1.0, 2.5);
  REQUIRE(a == Position(0.0, 0.0, 0.0));
}

TEST_CASE("Test linear time adjustment")
{
  // Typical use case
  double t = 3.0;
  Position a = Position(-1.0, -1.0, -1.0);
  Position b = Position(1.0, 1.0, 1.0);
  Position c = Position(0.0, 0.0, 0.0);
  _adjust_time(t, 2.0, a, b, c);
  REQUIRE(t == 2.5);
}
