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

} // namespace openmc
