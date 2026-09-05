#include <cmath>
#include <string>

#include <catch2/matchers/catch_matchers_templated.hpp>

#include "fmt/format.h"
#include "openmc/position.h"

using namespace openmc;

class EqualsPositionMatcher : public Catch::Matchers::MatcherGenericBase {
public:
  EqualsPositionMatcher(const Position& target, double margin = 1.0e-10)
    : target(target), margin(margin)
  {}

  bool match(const Position& other) const
  {
    return std::abs(other.x - target.x) <= margin &&
           std::abs(other.y - target.y) <= margin &&
           std::abs(other.z - target.z) <= margin;
  }

  std::string describe() const override
  {
    return fmt::format("Equals: ({}, {}, {})", target.x, target.y, target.z);
  }

private:
  const Position& target;
  double margin;
};

auto EqualsPosition(const Position& target, double margin)
  -> EqualsPositionMatcher
{
  return EqualsPositionMatcher {target, margin};
}
