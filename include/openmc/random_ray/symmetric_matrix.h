#ifndef OPENMC_SYMMETRIC_MATRIX_H
#define OPENMC_SYMMETRIC_MATRIX_H

#include <array>

#include "openmc/position.h"

namespace openmc {

// We make a typedef for the MomentArray type, which is a 3D vector of doubles.
using MomentArray = Position;

// The SymmetricMatrix class is a sparese representation a 3x3 symmetric
// matrix, with elements labeled as follows:
//
// | a b c |
// | b d e |
// | c e f |
//
class SymmetricMatrix {
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
  SymmetricMatrix() = default;
  SymmetricMatrix(double a, double b, double c, double d, double e, double f)
    : a(a), b(b), c(c), d(d), e(e), f(f)
  {}

  //----------------------------------------------------------------------------
  // Methods
  SymmetricMatrix inverse() const;
  double determinant() const;
  void set_to_zero() { a = b = c = d = e = f = 0; }
  void scale(double x)
  {
    a *= x;
    b *= x;
    c *= x;
    d *= x;
    e *= x;
    f *= x;
  }

  SymmetricMatrix& operator*= (double x)
  {
    a *= x;
    b *= x;
    c *= x;
    d *= x;
    e *= x;
    f *= x;
    return *this;
  }

  SymmetricMatrix& operator+=(const SymmetricMatrix& rhs)
  {
    a += rhs.a;
    b += rhs.b;
    c += rhs.c;
    d += rhs.d;
    e += rhs.e;
    f += rhs.f;
    return *this;
  }

  Position operator*(const Position& rhs) const
  {
    return {a * rhs.x + b * rhs.y + c * rhs.z,
            b * rhs.x + d * rhs.y + e * rhs.z,
            c * rhs.x + e * rhs.y + f * rhs.z};
  }

  std::array<double, 3> solve(const std::array<double, 3>& y) const;
  Position solve( const Position& y) const
  {
    std::array<double, 3> y_array = {y.x, y.y, y.z};
    std::array<double, 3> x_array = solve(y_array);
    return {x_array[0], x_array[1], x_array[2]};
  }

  void compute_spatial_moments_matrix(const Position& r, const Direction& u, const double& distance);
};

} // namespace openmc

#endif // OPENMC_SYMMETRIC_MATRIX_H