#include "openmc/random_ray/symmetric_matrix.h"

#include<cmath>

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
        return inv;
    }

    // Compute the adjoint matrix
    inv.a = d*f - e*e;
    inv.b = c*e - b*f;
    inv.c = b*e - c*d;
    inv.d = a*f - c*c;
    inv.e = b*c - a*e;
    inv.f = a*d - b*b;

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
    return a * (d*f - e*e) - b * (b*f - c*e) + c * (b*e - c*d);
}

} // namespace openmc