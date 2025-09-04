//! \file math_functions.h
//! A collection of elementary math functions.

#ifndef OPENMC_MATH_FUNCTIONS_H
#define OPENMC_MATH_FUNCTIONS_H

#include <cmath>
#include <complex>
#include <cstdlib>

#include "openmc/position.h"

namespace openmc {

//==============================================================================
//! Calculate the percentile of the standard normal distribution with a
//! specified probability level.
//!
//! \param p The probability level
//! \return The requested percentile
//==============================================================================

extern "C" double normal_percentile(double p);

//==============================================================================
//! Calculate the percentile of the Student's t distribution with a specified
//! probability level and number of degrees of freedom.
//!
//! \param p  The probability level
//! \param df The degrees of freedom
//! \return The requested percentile
//==============================================================================

extern "C" double t_percentile(double p, int df);

//==============================================================================
//! Calculate the n-th order Legendre polynomials at the value of x.
//!
//! \param n   The maximum order requested
//! \param x   The value to evaluate at; x is expected to be within [-1,1]
//! \param pnx The requested Legendre polynomials of order 0 to n (inclusive)
//!   evaluated at x.
//==============================================================================

extern "C" void calc_pn_c(int n, double x, double pnx[]);

//==============================================================================
//! Find the value of f(x) given a set of Legendre coefficients and the value
//! of x.
//!
//! \param n    The maximum order of the expansion
//! \param data The polynomial expansion coefficient data; without the (2l+1)/2
//!   factor.
//! \param x    The value to evaluate at; x is expected to be within [-1,1]
//! \return The requested Legendre polynomials of order 0 to n (inclusive)
//!   evaluated at x
//==============================================================================

extern "C" double evaluate_legendre(int n, const double data[], double x);

//==============================================================================
//! Calculate the n-th order real spherical harmonics for a given angle (in
//! terms of (u,v,w)) for all 0<=n and -m<=n<=n.
//!
//! \param n      The maximum order requested
//! \param uvw[3] The direction the harmonics are requested at
//! \param rn     The requested harmonics of order 0 to n (inclusive)
//!   evaluated at uvw.
//==============================================================================

extern "C" void calc_rn_c(int n, const double uvw[3], double rn[]);

void calc_rn(int n, Direction u, double rn[]);

//==============================================================================
//! Calculate the n-th order modified Zernike polynomial moment for a given
//! angle (rho, theta) location on the unit disk.
//!
//! This procedure uses the modified Kintner's method for calculating Zernike
//! polynomials as outlined in Chong, C. W., Raveendran, P., & Mukundan,
//! R. (2003). A comparative analysis of algorithms for fast computation of
//! Zernike moments. Pattern Recognition, 36(3), 731-742.
//! The normalization of the polynomials is such that the integral of Z_pq^2
//! over the unit disk is exactly pi.
//!
//! \param n   The maximum order requested
//! \param rho The radial parameter to specify location on the unit disk
//! \param phi The angle parameter to specify location on the unit disk
//! \param zn  The requested moments of order 0 to n (inclusive)
//!   evaluated at rho and phi.
//==============================================================================

extern "C" void calc_zn(int n, double rho, double phi, double zn[]);

//==============================================================================
//! Calculate only the even radial components of n-th order modified Zernike
//! polynomial moment with azimuthal dependency m = 0 for a given angle
//! (rho, theta) location on the unit disk.
//!
//! Since m = 0, n could only be even orders. Z_q0 = R_q0
//!
//! This procedure uses the modified Kintner's method for calculating Zernike
//! polynomials as outlined in Chong, C. W., Raveendran, P., & Mukundan,
//! R. (2003). A comparative analysis of algorithms for fast computation of
//! Zernike moments. Pattern Recognition, 36(3), 731-742.
//! The normalization of the polynomials is such that the integral of Z_pq^2
//! over the unit disk is exactly pi.
//!
//! \param n       The maximum order requested
//! \param rho     The radial parameter to specify location on the unit disk
//! \param phi     The angle parameter to specify location on the unit disk
//! \param zn_rad  The requested moments of order 0 to n (inclusive)
//!   evaluated at rho and phi when m = 0.
//==============================================================================

extern "C" void calc_zn_rad(int n, double rho, double zn_rad[]);

//==============================================================================
//! Rotate the direction cosines through a polar angle whose cosine is mu and
//! through an azimuthal angle sampled uniformly.
//!
//! This is done with direct sampling rather than rejection sampling as is done
//! in MCNP and Serpent.
//!
//! \param uvw[3] The initial, and final, direction vector
//! \param mu     The cosine of angle in lab or CM
//! \param phi    The azimuthal angle; will randomly chosen angle if a nullptr
//!   is passed
//! \param seed A pointer to the pseudorandom seed
//==============================================================================

extern "C" void rotate_angle_c(
  double uvw[3], double mu, const double* phi, uint64_t* seed);

Direction rotate_angle(
  Direction u, double mu, const double* phi, uint64_t* seed);

//==============================================================================
//! Constructs a natural cubic spline.
//!
//! Given a tabulated function y_i = f(x_i), this computes the second
//! derivative of the interpolating function at each x_i, which can then be
//! used in any subsequent calls to spline_interpolate or spline_integrate for
//! the same set of x and y values.
//!
//! \param n       Number of points
//! \param x       Values of the independent variable, which must be strictly
//!   increasing.
//! \param y       Values of the dependent variable.
//! \param[out] z  The second derivative of the interpolating function at each
//!   value of x.
//==============================================================================

void spline(int n, const double x[], const double y[], double z[]);

//==============================================================================
//! Determine the cubic spline interpolated y-value for a given x-value.
//!
//! \param n     Number of points
//! \param x     Values of the independent variable, which must be strictly
//!   increasing.
//! \param y     Values of the dependent variable.
//! \param z     The second derivative of the interpolating function at each
//!   value of x.
//! \param xint  Point at which to evaluate the cubic spline polynomial
//! \return      Interpolated value
//==============================================================================

double spline_interpolate(
  int n, const double x[], const double y[], const double z[], double xint);

//==============================================================================
//! Evaluate the definite integral of the interpolating cubic spline between
//! the given endpoints.
//!
//! \param n   Number of points
//! \param x   Values of the independent variable, which must be strictly
//!   increasing.
//! \param y   Values of the dependent variable.
//! \param z   The second derivative of the interpolating function at each
//!   value of x.
//! \param xa  Lower limit of integration
//! \param xb  Upper limit of integration
//! \return    Integral
//==============================================================================

double spline_integrate(int n, const double x[], const double y[],
  const double z[], double xa, double xb);

//! Evaluate the Faddeeva function
//!
//! \param z Complex argument
//! \return Faddeeva function evaluated at z
std::complex<double> faddeeva(std::complex<double> z);

//! Evaluate derivative of the Faddeeva function
//!
//! \param z Complex argument
//! \param order Order of the derivative
//! \return Derivative of Faddeeva function evaluated at z
std::complex<double> w_derivative(std::complex<double> z, int order);

//! Evaluate relative exponential function
//! \param x Real argument
//! \return (exp(x)-1)/x without loss of precision near 0
double exprel(double x);

//! Evaluate relative exponential function
//! \param x Real argument
//! \return log(1+x)/x without loss of precision near 0
double log1prel(double x);

//! Evaluate monomial difference
//! \param x Real argument
//! \param y Real argument
//! \param n Real argument
//! \return (x^n-y^n)/n without loss of precision near n=0
double monodiff(double x, double y, double n);

} // namespace openmc
#endif // OPENMC_MATH_FUNCTIONS_H
