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
  double xmin = -INFTY;
  double xmax = INFTY;
  double ymin = -INFTY;
  double ymax = INFTY;
  double zmin = -INFTY;
  double zmax = INFTY;

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
    xmin = std::max(xmin, other.xmin);
    xmax = std::min(xmax, other.xmax);
    ymin = std::max(ymin, other.ymin);
    ymax = std::min(ymax, other.ymax);
    zmin = std::max(zmin, other.zmin);
    zmax = std::min(zmax, other.zmax);
    return *this;
  }

  // union operator
  inline BoundingBox& operator|=(const BoundingBox& other)
  {
    xmin = std::min(xmin, other.xmin);
    xmax = std::max(xmax, other.xmax);
    ymin = std::min(ymin, other.ymin);
    ymax = std::max(ymax, other.ymax);
    zmin = std::min(zmin, other.zmin);
    zmax = std::max(zmax, other.zmax);
    return *this;
  }

  inline Position min() const { return {xmin, ymin, zmin}; }
  inline Position max() const { return {xmax, ymax, zmax}; }
};

} // namespace openmc

#endif
