#ifndef OPENMC_MOMENT_MATRIX_H
#define OPENMC_MOMENT_MATRIX_H

#include <array>

#include "openmc/position.h"

namespace openmc {

// The MomentArray class is a 3-element array representing the x, y, and z
// moments. It is defined as an alias for the Position class to allow for
// dot products and other operations with Position and Direction objects.
// TODO: This class could in theory have 32-bit instead of 64-bit FP values.
using MomentArray = Position;

// The MomentMatrix class is a sparese representation a 3x3 symmetric
// matrix, with elements labeled as follows:
//
// | a b c |
// | b d e |
// | c e f |
//
// This class uses FP64 values as objects will be accumulated to over many
// iterations.
class MomentMatrix {
public:
  //----------------------------------------------------------------------------
  // Public data members
  double a;
  double b;
  double c;
  double d;
  double e;
  double f;

  //----------------------------------------------------------------------------
  // Constructors
  MomentMatrix() = default;
  MomentMatrix(double a, double b, double c, double d, double e, double f)
    : a {a}, b {b}, c {c}, d {d}, e {e}, f {f}
  {}

  //----------------------------------------------------------------------------
  // Methods
  MomentMatrix inverse() const;
  double determinant() const;
  void set_to_zero() { a = b = c = d = e = f = 0; }

  inline MomentMatrix& operator*=(double x)
  {
    a *= x;
    b *= x;
    c *= x;
    d *= x;
    e *= x;
    f *= x;
    return *this;
  }

  inline MomentMatrix operator*(double x) const
  {
    MomentMatrix m_copy = *this;
    m_copy *= x;
    return m_copy;
  }

  inline MomentMatrix& operator+=(const MomentMatrix& rhs)
  {
    a += rhs.a;
    b += rhs.b;
    c += rhs.c;
    d += rhs.d;
    e += rhs.e;
    f += rhs.f;
    return *this;
  }

  MomentArray operator*(const MomentArray& rhs) const
  {
    return {a * rhs.x + b * rhs.y + c * rhs.z,
      b * rhs.x + d * rhs.y + e * rhs.z, c * rhs.x + e * rhs.y + f * rhs.z};
  }

  MomentArray solve(const MomentArray& y) const;

  void compute_spatial_moments_matrix(
    const Position& r, const Direction& u, const double& distance);
};

} // namespace openmc

#endif // OPENMC_MOMENT_MATRIX_H