#ifndef OPENMC_MATH_H
#define OPENMC_MATH_H

//! \file math.h
//! \brief Transcendental functions used by the transport code.
//!
//! Transport code calls transcendental functions through the openmc namespace
//! (openmc::exp, openmc::log, ...) rather than through std::. By default these
//! are the standard library functions, so nothing changes. When OpenMC is
//! configured with OPENMC_USE_COREMATH=ON they map instead to the correctly
//! rounded implementations of the CORE-MATH project, which return the same
//! bits on every platform and compiler. System math libraries (glibc, Apple
//! libSystem, the MSVC runtime) round these functions differently, and a
//! one-bit difference in a sampled collision distance is enough to send two
//! particle histories down different paths, which is why regression results
//! differ between operating systems.
//!
//! Routing the calls through one header is what makes the backend
//! selectable; it costs nothing when the standard library is used.

#include <cmath>

#ifdef OPENMC_COREMATH

extern "C" {
double cr_acos(double);
double cr_asin(double);
double cr_atan(double);
double cr_atan2(double, double);
double cr_cbrt(double);
double cr_cos(double);
double cr_cosh(double);
double cr_erf(double);
double cr_erfc(double);
double cr_exp(double);
double cr_expm1(double);
double cr_hypot(double, double);
double cr_lgamma(double);
double cr_log(double);
double cr_log1p(double);
double cr_log2(double);
double cr_pow(double, double);
double cr_sin(double);
double cr_sinh(double);
double cr_tan(double);
double cr_tgamma(double);
}

namespace openmc {

//! Arc cosine of x, correctly rounded
inline double acos(double x)
{
  return cr_acos(x);
}
//! Arc sine of x, correctly rounded
inline double asin(double x)
{
  return cr_asin(x);
}
//! Arc tangent of x, correctly rounded
inline double atan(double x)
{
  return cr_atan(x);
}
//! Angle of the point (x, y), in (-pi, pi], correctly rounded
inline double atan2(double y, double x)
{
  return cr_atan2(y, x);
}
//! Cube root of x, correctly rounded
inline double cbrt(double x)
{
  return cr_cbrt(x);
}
//! Cosine of x, correctly rounded
inline double cos(double x)
{
  return cr_cos(x);
}
//! Hyperbolic cosine of x, correctly rounded
inline double cosh(double x)
{
  return cr_cosh(x);
}
//! Error function of x, correctly rounded
inline double erf(double x)
{
  return cr_erf(x);
}
//! Complementary error function of x, correctly rounded
inline double erfc(double x)
{
  return cr_erfc(x);
}
//! Exponential of x, correctly rounded
inline double exp(double x)
{
  return cr_exp(x);
}
//! exp(x) - 1, accurate for small x of x, correctly rounded
inline double expm1(double x)
{
  return cr_expm1(x);
}
//! sqrt(x^2 + y^2) without intermediate overflow, correctly rounded
inline double hypot(double x, double y)
{
  return cr_hypot(x, y);
}
//! Natural logarithm of the absolute value of the gamma function of x,
//! correctly rounded
inline double lgamma(double x)
{
  return cr_lgamma(x);
}
//! Natural logarithm of x, correctly rounded
inline double log(double x)
{
  return cr_log(x);
}
//! log(1 + x), accurate for small x of x, correctly rounded
inline double log1p(double x)
{
  return cr_log1p(x);
}
//! Base-2 logarithm of x, correctly rounded
inline double log2(double x)
{
  return cr_log2(x);
}
//! x raised to the power y, correctly rounded
inline double pow(double x, double y)
{
  return cr_pow(x, y);
}
//! Sine of x, correctly rounded
inline double sin(double x)
{
  return cr_sin(x);
}
//! Hyperbolic sine of x, correctly rounded
inline double sinh(double x)
{
  return cr_sinh(x);
}
//! Tangent of x, correctly rounded
inline double tan(double x)
{
  return cr_tan(x);
}
//! Gamma function of x, correctly rounded
inline double tgamma(double x)
{
  return cr_tgamma(x);
}

} // namespace openmc

#else

namespace openmc {
using std::acos;
using std::asin;
using std::atan;
using std::atan2;
using std::cbrt;
using std::cos;
using std::cosh;
using std::erf;
using std::erfc;
using std::exp;
using std::expm1;
using std::hypot;
using std::lgamma;
using std::log;
using std::log1p;
using std::log2;
using std::pow;
using std::sin;
using std::sinh;
using std::tan;
using std::tgamma;
} // namespace openmc

#endif // OPENMC_COREMATH

#endif // OPENMC_MATH_H
