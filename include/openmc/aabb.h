#ifndef OPENMC_AABB_H
#define OPENMC_AABB_H

#include "position.h"

namespace openmc {

// This axis-aligned bounding box class (AABB) is designed to work easier with
// partitioners than openmc::BoundingBox
struct AABB {
  AABB();
  AABB(const Position& mi, const Position& ma);

  void extend_max(const Position& val);
  void extend_min(const Position& val);

  void extend(const Position& pos);
  void extend(const AABB& other_box);

  Position get_center() const;
  double surface_area() const;
  double volume() const;
  bool contains(const Position& pos) const;

  void reset();

  bool operator==(const AABB& other) const;
  bool operator!=(const AABB& other) const;

  Position min_;
  Position max_;
};

}; // namespace openmc

#endif