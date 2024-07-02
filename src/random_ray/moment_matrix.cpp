#include "openmc/error.h"
#include "openmc/random_ray/moment_matrix.h"

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
MomentMatrix MomentMatrix::inverse() const
{
  MomentMatrix inv;

  // Check if the determinant is zero
  double det = determinant();
  if (det < std::abs(1.0e-10)) {
    // Set the inverse to zero. In effect, this will
    // result in all the linear terms of the source becoming
    // zero, leaving just the flat source. 
    inv.set_to_zero();
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
  inv *= 1.0 / det;

  return inv;
}

// Computes the determinant of a 3x3 symmetric
// matrix, with elements labeled as follows:
//
// | a b c |
// | b d e |
// | c e f |
double MomentMatrix::determinant() const
{
  return a * (d * f - e * e) - b * (b * f - c * e) + c * (b * e - c * d);
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
void MomentMatrix::compute_spatial_moments_matrix(
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