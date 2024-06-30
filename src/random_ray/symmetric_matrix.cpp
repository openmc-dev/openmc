#include "openmc/random_ray/symmetric_matrix.h"
#include "openmc/error.h"

#include <cmath>

namespace openmc {

//==============================================================================
// UpperTriangular implementation
//==============================================================================

// Inters the following 3x3 smmetric matrix labeled as:
//
// | a b c |
// | b d e |
// | c e f |
//
// We first check the determinant to ensure it is non-zero
// before proceeding with the inversion. If the determinant
// is zero, we return a matrix of zeros.
// Inversion is calculated by computing the adjoint matrix
// first, and then the inverse can be computed as:
// A^-1  = 1/det(A) * adj(A)
SymmetricMatrix SymmetricMatrix::inverse() const
{
  SymmetricMatrix inv;

  // Check if the determinant is zero
  double det = determinant();
  if (det < std::abs(1.0e-10)) {
    inv.set_to_zero();
    // fatal_error("Matrix is singular and cannot be inverted.");
    return inv;
  }

  // Compute the adjoint matrix
  inv.a = d * f - e * e;
  inv.b = c * e - b * f;
  inv.c = b * e - c * d;
  inv.d = a * f - c * c;
  inv.e = b * c - a * e;
  inv.f = a * d - b * b;

  // A^-1 = 1/det(A) * adj(A)
  inv.scale(1.0 / det);

  return inv;
}

// Computes the determinant of a 3x3 symmetric
// matrix, with elements labeled as follows:
//
// | a b c |
// | b d e |
// | c e f |
double SymmetricMatrix::determinant() const
{
  return a * (d * f - e * e) - b * (b * f - c * e) + c * (b * e - c * d);
}

std::array<double, 3> SymmetricMatrix::solve(
  const std::array<double, 3>& y) const
{
  // Check for positive definiteness and calculate Cholesky decomposition
  if (a <= 1e-10) {
    fatal_error("Matrix is not positive definite (element 'a' non-positive).");
  }

  double L11 = std::sqrt(a);
  double L21 = b / L11; // Using matrix element 'b' directly
  double L31 = c / L11;

  double tmp = d - L21 * L21;
  if (tmp <= 1e-10) {
    fatal_error("Matrix is not positive definite at second diagonal element.");
  }
  double L22 = std::sqrt(tmp);

  double L32 = (e - L21 * L31) / L22;

  tmp = f - L31 * L31 - L32 * L32;
  if (tmp <= 1e-10) {
    fatal_error("Matrix is not positive definite at third diagonal element.");
  }
  double L33 = std::sqrt(tmp);

  // Solving Ly = y (forward substitution)
  double y1 = y[0] / L11;
  double y2 = (y[1] - L21 * y1) / L22;
  double y3 = (y[2] - L31 * y1 - L32 * y2) / L33;

  // Solving L^T x = y (backward substitution)
  double x3 = y3 / L33;
  double x2 = (y2 - L32 * x3) / L22;
  double x1 = (y1 - L21 * x2 - L31 * x3) / L11;

  return {x1, x2, x3};
}

// Compute a 3x3 spatial moment matrix based on a single ray crossing.
// The matrix is symmetric, and is defined as:
//
// | a b c |
// | b d e |
// | c e f |
//
// The estimate of the obect's spatial moments matrix is computed based on the
// midpoint of the ray's crossing, the direction of the ray, and the distance
// the ray traveled through the 3D object.
void SymmetricMatrix::compute_spatial_moments_matrix(
  const Position& r, const Direction& u, const double& distance)
{
  constexpr double one_over_twelve = 1.0 / 12.0;
  const double distance2_12 = distance * distance * one_over_twelve;
  a = r[0] * r[0] + u[0] * u[0] * distance2_12;
  b = r[0] * r[1] + u[0] * u[1] * distance2_12;
  c = r[0] * r[2] + u[0] * u[2] * distance2_12;
  d = r[1] * r[1] + u[1] * u[1] * distance2_12;
  e = r[1] * r[2] + u[1] * u[2] * distance2_12;
  f = r[2] * r[2] + u[2] * u[2] * distance2_12;
}

} // namespace openmc