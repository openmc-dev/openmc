#ifndef OPENMC_BOUNDING_BOX_H
#define OPENMC_BOUNDING_BOX_H

#include "openmc/constants.h"
#include "openmc/position.h"

namespace openmc {

//==============================================================================
//! Coordinates for an axis-aligned cuboid bounds a geometric object.
//==============================================================================

struct BoundingBox
{
  double xmin = -INFTY;
  double xmax = INFTY;
  double ymin = -INFTY;
  double ymax = INFTY;
  double zmin = -INFTY;
  double zmax = INFTY;


  inline BoundingBox operator &(const BoundingBox& other) {
    BoundingBox result = *this;
    return result &= other;
  }

  inline BoundingBox operator |(const BoundingBox& other) {
    BoundingBox result = *this;
    return result |= other;
  }

  // intersect operator
  inline BoundingBox& operator &=(const BoundingBox& other) {
    xmin = std::max(xmin, other.xmin);
    xmax = std::min(xmax, other.xmax);
    ymin = std::max(ymin, other.ymin);
    ymax = std::min(ymax, other.ymax);
    zmin = std::max(zmin, other.zmin);
    zmax = std::min(zmax, other.zmax);
    return *this;
  }

  // union operator
  inline BoundingBox& operator |=(const BoundingBox& other) {
    xmin = std::min(xmin, other.xmin);
    xmax = std::max(xmax, other.xmax);
    ymin = std::min(ymin, other.ymin);
    ymax = std::max(ymax, other.ymax);
    zmin = std::min(zmin, other.zmin);
    zmax = std::max(zmax, other.zmax);
    return *this;
  }

  // ensure bounding box contains a point
  inline void update(const Position& r) {
    xmin = std::min(xmin, r.x);
    xmax = std::max(xmax, r.x);
    ymin = std::min(ymin, r.y);
    ymax = std::max(ymax, r.y);
    zmin = std::min(zmin, r.z);
    zmax = std::max(zmax, r.z);
  }

  // check if a point is in the box
  inline bool contains(const Position& r) const {
    if (r.x < xmin || r.x > xmax) { return false; }
    if (r.y < ymin || r.y > ymax) { return false; }
    if (r.z < zmin || r.z > zmax) { return false; }
    return true;
  }

};

} // namespace openmc
#endif // OPENMC_BOUNDING_BOX_H
