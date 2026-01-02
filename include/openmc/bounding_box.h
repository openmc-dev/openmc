#ifndef OPENMC_BOUNDING_BOX_H
#define OPENMC_BOUNDING_BOX_H

#include <algorithm> // for min, max

#include "openmc/constants.h"
#include "openmc/position.h"

namespace openmc {

//==============================================================================
//! Coordinates for an axis-aligned cuboid that bounds a geometric object.
//==============================================================================

struct BoundingBox {
  Position min = {-INFTY, -INFTY, -INFTY};
  Position max = {INFTY, INFTY, INFTY};

  // Constructors
  BoundingBox() = default;
  BoundingBox(Position min_, Position max_) : min {min_}, max {max_} {}

  // Static factory methods
  static BoundingBox infinite() { return {}; }
  static BoundingBox inverted()
  {
    return {{INFTY, INFTY, INFTY}, {-INFTY, -INFTY, -INFTY}};
  }

  inline BoundingBox operator&(const BoundingBox& other)
  {
    BoundingBox result = *this;
    return result &= other;
  }

  inline BoundingBox operator|(const BoundingBox& other)
  {
    BoundingBox result = *this;
    return result |= other;
  }

  // intersect operator
  inline BoundingBox& operator&=(const BoundingBox& other)
  {
    min.x = std::max(min.x, other.min.x);
    min.y = std::max(min.y, other.min.y);
    min.z = std::max(min.z, other.min.z);
    max.x = std::min(max.x, other.max.x);
    max.y = std::min(max.y, other.max.y);
    max.z = std::min(max.z, other.max.z);
    return *this;
  }

  // union operator
  inline BoundingBox& operator|=(const BoundingBox& other)
  {
    min.x = std::min(min.x, other.min.x);
    min.y = std::min(min.y, other.min.y);
    min.z = std::min(min.z, other.min.z);
    max.x = std::max(max.x, other.max.x);
    max.y = std::max(max.y, other.max.y);
    max.z = std::max(max.z, other.max.z);
    return *this;
  }
};

} // namespace openmc

#endif
