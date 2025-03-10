#ifndef OPENMC_POSITION_H
#define OPENMC_POSITION_H

#include <cmath> // for sqrt
#include <iostream>
#include <stdexcept> // for out_of_range

#include "fmt/format.h"
#include "openmc/array.h"
#include "openmc/vector.h"

namespace openmc {

//==============================================================================
//! Type representing a position in Cartesian coordinates
//==============================================================================

struct Position {
  // Constructors
  Position() = default;
  Position(double x_, double y_, double z_) : x {x_}, y {y_}, z {z_} {};
  Position(const double xyz[]) : x {xyz[0]}, y {xyz[1]}, z {xyz[2]} {};
  Position(const vector<double>& xyz) : x {xyz[0]}, y {xyz[1]}, z {xyz[2]} {};
  Position(const array<double, 3>& xyz) : x {xyz[0]}, y {xyz[1]}, z {xyz[2]} {};

  // Unary operators
  Position& operator+=(Position);
  Position& operator+=(double);
  Position& operator-=(Position);
  Position& operator-=(double);
  Position& operator*=(Position);
  Position& operator*=(double);
  Position& operator/=(Position);
  Position& operator/=(double);
  Position operator-() const;

  const double& operator[](int i) const
  {
    switch (i) {
    case 0:
      return x;
    case 1:
      return y;
    case 2:
      return z;
    default:
      throw std::out_of_range {"Index in Position must be between 0 and 2."};
    }
  }
  double& operator[](int i)
  {
    switch (i) {
    case 0:
      return x;
    case 1:
      return y;
    case 2:
      return z;
    default:
      throw std::out_of_range {"Index in Position must be between 0 and 2."};
    }
  }

  // Access to x, y, or z by compile time known index (specializations below)
  template<int i>
  const double& get() const
  {
    throw std::out_of_range {"Index in Position must be between 0 and 2."};
  }
  template<int i>
  double& get()
  {
    throw std::out_of_range {"Index in Position must be between 0 and 2."};
  }

  // Other member functions

  //! Dot product of two vectors
  //! \param[in] other Vector to take dot product with
  //! \result Resulting dot product
  inline double dot(Position other) const
  {
    return x * other.x + y * other.y + z * other.z;
  }
  inline double norm() const { return std::sqrt(x * x + y * y + z * z); }
  inline Position cross(Position other) const
  {
    return {y * other.z - z * other.y, z * other.x - x * other.z,
      x * other.y - y * other.x};
  }

  //! Reflect a direction across a normal vector
  //! \param[in] other Vector to reflect across
  //! \result Reflected vector
  Position reflect(Position n) const;

  //! Rotate the position by applying a rotation matrix
  template<typename T>
  Position rotate(const T& rotation) const
  {
    return {x * rotation[0] + y * rotation[1] + z * rotation[2],
      x * rotation[3] + y * rotation[4] + z * rotation[5],
      x * rotation[6] + y * rotation[7] + z * rotation[8]};
  }

  //! Rotate the position by applying the inverse of a rotation matrix
  //! using the fact that rotation matrices are orthonormal.
  template<typename T>
  Position inverse_rotate(const T& rotation) const
  {
    return {x * rotation[0] + y * rotation[3] + z * rotation[6],
      x * rotation[1] + y * rotation[4] + z * rotation[7],
      x * rotation[2] + y * rotation[5] + z * rotation[8]};
  }

  // Data members
  double x = 0.;
  double y = 0.;
  double z = 0.;
};

// Compile-time known member index access functions
template<>
inline const double& Position::get<0>() const
{
  return x;
}
template<>
inline const double& Position::get<1>() const
{
  return y;
}
template<>
inline const double& Position::get<2>() const
{
  return z;
}
template<>
inline double& Position::get<0>()
{
  return x;
}
template<>
inline double& Position::get<1>()
{
  return y;
}
template<>
inline double& Position::get<2>()
{
  return z;
}

// Binary operators
inline Position operator+(Position a, Position b)
{
  return a += b;
}
inline Position operator+(Position a, double b)
{
  return a += b;
}
inline Position operator+(double a, Position b)
{
  return b += a;
}

inline Position operator-(Position a, Position b)
{
  return a -= b;
}
inline Position operator-(Position a, double b)
{
  return a -= b;
}
inline Position operator-(double a, Position b)
{
  return b -= a;
}

inline Position operator*(Position a, Position b)
{
  return a *= b;
}
inline Position operator*(Position a, double b)
{
  return a *= b;
}
inline Position operator*(double a, Position b)
{
  return b *= a;
}

inline Position operator/(Position a, Position b)
{
  return a /= b;
}
inline Position operator/(Position a, double b)
{
  return a /= b;
}
inline Position operator/(double a, Position b)
{
  return b /= a;
}

inline Position Position::reflect(Position n) const
{
  const double projection = n.dot(*this);
  const double magnitude = n.dot(n);
  n *= (2.0 * projection / magnitude);
  return *this - n;
}

inline bool operator==(Position a, Position b)
{
  return a.x == b.x && a.y == b.y && a.z == b.z;
}

inline bool operator!=(Position a, Position b)
{
  return a.x != b.x || a.y != b.y || a.z != b.z;
}

std::ostream& operator<<(std::ostream& os, Position a);

//==============================================================================
//! Type representing a vector direction in Cartesian coordinates
//==============================================================================

using Direction = Position;

} // namespace openmc

namespace fmt {

template<>
struct formatter<openmc::Position> : formatter<std::string> {
  template<typename FormatContext>
#if FMT_VERSION >= 110000 // Version 11.0.0 and above
  auto format(const openmc::Position& pos, FormatContext& ctx) const {
#else // For versions below 11.0.0
  auto format(const openmc::Position& pos, FormatContext& ctx)
  {
#endif
    return formatter<std::string>::format(
      fmt::format("({}, {}, {})", pos.x, pos.y, pos.z), ctx);
}
}; // namespace fmt

} // namespace fmt

#endif // OPENMC_POSITION_H
