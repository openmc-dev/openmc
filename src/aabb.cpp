#include "openmc/aabb.h"
#include <float.h>

namespace openmc {

// Take the minimum of two vectors
inline Position vmin(Position lhs, Position rhs)
{
  Position res;
  res.x = std::min(lhs.x, rhs.x);
  res.y = std::min(lhs.y, rhs.y);
  res.z = std::min(lhs.z, rhs.z);
  return res;
}

// Take the maximum of two vectors
inline Position vmax(Position lhs, Position rhs)
{
  Position res;
  res.x = std::max(lhs.x, rhs.x);
  res.y = std::max(lhs.y, rhs.y);
  res.z = std::max(lhs.z, rhs.z);
  return res;
}

// Initialize
AABB::AABB()
  : min_(DBL_MAX, DBL_MAX, DBL_MAX), max_(-DBL_MAX, -DBL_MAX, -DBL_MAX)
{}

// Initialize
AABB::AABB(const Position& mi, const Position& ma) : min_(mi), max_(ma) {}

// Set everything to infinity
void AABB::reset()
{
  AABB clear;
  *this = clear;
}

void AABB::extend_max(const Position& val)
{
  max_ = vmax(max_, val);
}

void AABB::extend_min(const Position& val)
{
  min_ = vmin(min_, val);
}

void AABB::extend(const Position& pos)
{
  extend_max(pos);
  extend_min(pos);
}

double AABB::surface_area() const
{
  Position side_lengths = max_ - min_;

  return 2 * (side_lengths.x * (side_lengths.y + side_lengths.z) +
               side_lengths.y * side_lengths.z);
}

double AABB::volume() const
{
  Position side_lengths = max_ - min_;
  return side_lengths.x * side_lengths.y * side_lengths.z;
}

void AABB::extend(const AABB& other_box)
{
  extend_max(other_box.max_);
  extend_min(other_box.min_);
}

Position AABB::get_center() const
{
  return (min_ + max_) * 0.5f;
}

bool AABB::contains(const Position& pos) const
{
  return (min_.x <= pos.x && min_.y <= pos.y && min_.z <= pos.z &&
          max_.x >= pos.x && max_.y >= pos.y && max_.z >= pos.z);
}

bool AABB::operator==(const AABB& other) const
{
  return (min_.x == other.min_.x && min_.y == other.min_.y &&
          min_.z == other.min_.z && max_.x == other.max_.x &&
          max_.y == other.max_.y && max_.z == other.max_.z);
}

bool AABB::operator!=(const AABB& other) const
{
  return !(*this == other);
}

}; // namespace openmc