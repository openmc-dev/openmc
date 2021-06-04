#ifndef OPENMC_POSITION_H
#define OPENMC_POSITION_H

#include <cmath> // for sqrt
#include <iostream>
#include <stdexcept> // for out_of_range

#include "openmc/array.h"
#include "openmc/vector.h"

namespace openmc {

//==============================================================================
//! Type representing a position in Cartesian coordinates
//==============================================================================

struct Position {
  // Constructors
  Position() = default;
  HD Position(double x_, double y_, double z_) : x {x_}, y {y_}, z {z_} {};
  HD Position(const double xyz[]) : x {xyz[0]}, y {xyz[1]}, z {xyz[2]} {};
  HD Position(const vector<double>& xyz)
    : x {xyz[0]}, y {xyz[1]}, z {xyz[2]} {};
  Position(const vector<float>& xyz)
    : x {xyz[0]}, y {xyz[1]}, z {xyz[2]} {};
  Position(const array<double, 3>& xyz) : x{xyz[0]}, y{xyz[1]}, z{xyz[2]} { };
  Position(const std::array<double, 3>& xyz)
    : x {xyz[0]}, y {xyz[1]}, z {xyz[2]} {};

  // Unary operators
  HD Position& operator+=(Position);
  HD Position& operator+=(double);
  HD Position& operator-=(Position);
  HD Position& operator-=(double);
  HD Position& operator*=(Position);
  HD Position& operator*=(double);
  HD Position& operator/=(Position);
  HD Position& operator/=(double);
  HD Position operator-() const;

  HD const double& operator[](int i) const
  {
    switch (i) {
      case 0: return x;
      case 1: return y;
      case 2: return z;
#ifndef __CUDA_ARCH__
      default:
        throw std::out_of_range{"Index in Position must be between 0 and 2."};
#else
      default:
        asm("trap;");
        return x; // suppress warning
#endif
      }
  }
  HD double& operator[](int i)
  {
    switch (i) {
      case 0: return x;
      case 1: return y;
      case 2: return z;
#ifndef __CUDA_ARCH__
      default:
        throw std::out_of_range{"Index in Position must be between 0 and 2."};
#else
      default:
        asm("trap;");
        return x; // suppress warning
#endif
      }
  }

  // Access to x, y, or z by compile time known index (specializations below)
  template<int i>
  HD const double& get() const
  {
    throw std::out_of_range {"Index in Position must be between 0 and 2."};
  }
  template<int i>
  HD double& get()
  {
    throw std::out_of_range {"Index in Position must be between 0 and 2."};
  }

  // Other member functions

  //! Dot product of two vectors
  //! \param[in] other Vector to take dot product with
  //! \result Resulting dot product
  HD inline double dot(Position other) const
  {
    return x*other.x + y*other.y + z*other.z;
  }
  HD inline double norm() const { return std::sqrt(x * x + y * y + z * z); }

  //! Reflect a direction across a normal vector
  //! \param[in] other Vector to reflect across
  //! \result Reflected vector
  HD Position reflect(Position n) const;

  //! Rotate the position based on a rotation matrix
  HD Position rotate(const vector<double>& rotation) const;

  // Data members
  double x = 0.;
  double y = 0.;
  double z = 0.;
};

// Compile-time known member index access functions
template<>
HD inline const double& Position::get<0>() const
{
  return x;
}
template<>
HD inline const double& Position::get<1>() const
{
  return y;
}
template<>
HD inline const double& Position::get<2>() const
{
  return z;
}
template<>
HD inline double& Position::get<0>()
{
  return x;
}
template<>
HD inline double& Position::get<1>()
{
  return y;
}
template<>
HD inline double& Position::get<2>()
{
  return z;
}

// Binary operators
HD inline Position operator+(Position a, Position b)
{
  return a += b;
}
HD inline Position operator+(Position a, double b)
{
  return a += b;
}
HD inline Position operator+(double a, Position b)
{
  return b += a;
}

HD inline Position operator-(Position a, Position b)
{
  return a -= b;
}
HD inline Position operator-(Position a, double b)
{
  return a -= b;
}
HD inline Position operator-(double a, Position b)
{
  return b -= a;
}

HD inline Position operator*(Position a, Position b)
{
  return a *= b;
}
HD inline Position operator*(Position a, double b)
{
  return a *= b;
}
HD inline Position operator*(double a, Position b)
{
  return b *= a;
}

HD inline Position operator/(Position a, Position b)
{
  return a /= b;
}
HD inline Position operator/(Position a, double b)
{
  return a /= b;
}
HD inline Position operator/(double a, Position b)
{
  return b /= a;
}

HD inline Position Position::reflect(Position n) const
{
  const double projection = n.dot(*this);
  const double magnitude = n.dot(n);
  n *= (2.0 * projection / magnitude);
  return *this - n;
}

HD inline bool operator==(Position a, Position b)
{return a.x == b.x && a.y == b.y && a.z == b.z;}

HD inline bool operator!=(Position a, Position b)
{return a.x != b.x || a.y != b.y || a.z != b.z;}

std::ostream& operator<<(std::ostream& os, Position a);

//==============================================================================
//! Type representing a vector direction in Cartesian coordinates
//==============================================================================

using Direction = Position;

} // namespace openmc

#endif // OPENMC_POSITION_H
