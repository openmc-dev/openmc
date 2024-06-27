#ifndef OPENMC_SYMMETRIC_MATRIX_H
#define OPENMC_SYMMETRIC_MATRIX_H

namespace openmc {

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
};

} // namespace openmc

#endif // OPENMC_SYMMETRIC_MATRIX_H