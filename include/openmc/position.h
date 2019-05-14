#ifndef OPENMC_POSITION_H
#define OPENMC_POSITION_H

#include <array>
#include <cmath> // for sqrt
#include <iostream>
#include <stdexcept> // for out_of_range
#include <vector>

namespace openmc {

//==============================================================================
//! Type representing a position in Cartesian coordinates
//==============================================================================

struct Position {
  // Constructors
  Position() = default;
  Position(double x_, double y_, double z_) : x{x_}, y{y_}, z{z_} { };
  Position(const double xyz[]) : x{xyz[0]}, y{xyz[1]}, z{xyz[2]} { };
  Position(const std::vector<double>& xyz) : x{xyz[0]}, y{xyz[1]}, z{xyz[2]} { };
  Position(const std::array<double, 3>& xyz) : x{xyz[0]}, y{xyz[1]}, z{xyz[2]} { };

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

  const double& operator[](int i) const {
    switch (i) {
      case 0: return x;
      case 1: return y;
      case 2: return z;
      default:
        throw std::out_of_range{"Index in Position must be between 0 and 2."};
    }
  }
  double& operator[](int i) {
    switch (i) {
      case 0: return x;
      case 1: return y;
      case 2: return z;
      default:
        throw std::out_of_range{"Index in Position must be between 0 and 2."};
    }
  }

  // Other member functions

  //! Dot product of two vectors
  //! \param[in] other Vector to take dot product with
  //! \result Resulting dot product
  inline double dot(Position other) const {
    return x*other.x + y*other.y + z*other.z;
  }
  inline double norm() const {
    return std::sqrt(x*x + y*y + z*z);
  }

  //! Rotate the position based on a rotation matrix
  Position rotate(const std::vector<double>& rotation) const;

  // Data members
  double x = 0.;
  double y = 0.;
  double z = 0.;
};

// Binary operators
inline Position operator+(Position a, Position b) { return a += b; }
inline Position operator+(Position a, double b)   { return a += b; }
inline Position operator+(double a, Position b)   { return b += a; }

inline Position operator-(Position a, Position b) { return a -= b; }
inline Position operator-(Position a, double b)   { return a -= b; }
inline Position operator-(double a, Position b)   { return b -= a; }

inline Position operator*(Position a, Position b) { return a *= b; }
inline Position operator*(Position a, double b)   { return a *= b; }
inline Position operator*(double a, Position b)   { return b *= a; }

inline Position operator/(Position a, Position b) { return a /= b; }
inline Position operator/(Position a, double b)   { return a /= b; }
inline Position operator/(double a, Position b)   { return b /= a; }

inline bool operator==(Position a, Position b)
{return a.x == b.x && a.y == b.y && a.z == b.z;}

inline bool operator!=(Position a, Position b)
{return a.x != b.x || a.y != b.y || a.z != b.z;}

std::ostream& operator<<(std::ostream& os, Position a);

//==============================================================================
//! Type representing a vector direction in Cartesian coordinates
//==============================================================================

using Direction = Position;

} // namespace openmc

#endif // OPENMC_POSITION_H
