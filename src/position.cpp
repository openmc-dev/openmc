#include "openmc/position.h"

namespace openmc {

//==============================================================================
// Position implementation
//==============================================================================

Position&
Position::operator+=(Position other)
{
  x += other.x;
  y += other.y;
  z += other.z;
  return *this;
}

Position&
Position::operator+=(double v)
{
  x += v;
  y += v;
  z += v;
  return *this;
}

Position&
Position::operator-=(Position other)
{
  x -= other.x;
  y -= other.y;
  z -= other.z;
  return *this;
}

Position&
Position::operator-=(double v)
{
  x -= v;
  y -= v;
  z -= v;
  return *this;
}

Position&
Position::operator*=(Position other)
{
  x *= other.x;
  y *= other.y;
  z *= other.z;
  return *this;
}

Position&
Position::operator*=(double v)
{
  x *= v;
  y *= v;
  z *= v;
  return *this;
}

Position&
Position::operator/=(Position other)
{
  x /= other.x;
  y /= other.y;
  z /= other.z;
  return *this;
}

Position&
Position::operator/=(double v)
{
  x /= v;
  y /= v;
  z /= v;
  return *this;
}

Position
Position::operator-() const
{
  return {-x, -y, -z};
}

Position
Position::rotate(const std::vector<double>& rotation) const
{
  return {
    x*rotation[0] + y*rotation[1] + z*rotation[2],
    x*rotation[3] + y*rotation[4] + z*rotation[5],
    x*rotation[6] + y*rotation[7] + z*rotation[8]
  };
}

std::ostream&
operator<<(std::ostream& os, Position r)
{
  os << "(" << r.x << ", " << r.y << ", " << r.z << ")";
  return os;
}

} // namespace openmc
