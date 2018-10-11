//! \file math_functions.h
//! A collection of elementary math functions.

#ifndef OPENMC_MATH_FUNCTIONS_H
#define OPENMC_MATH_FUNCTIONS_H

#include <cmath>
#include <cstdlib>

#include "openmc/constants.h"
#include "openmc/position.h"
#include "openmc/random_lcg.h"


namespace openmc {

//==============================================================================
//! Calculate the percentile of the standard normal distribution with a
//! specified probability level.
//!
//! @param p The probability level
//! @return The requested percentile
//==============================================================================

extern "C" double normal_percentile_c(double p);

//==============================================================================
//! Calculate the percentile of the Student's t distribution with a specified
//! probability level and number of degrees of freedom.
//!
//! @param p  The probability level
//! @param df The degrees of freedom
//! @return The requested percentile
//==============================================================================

extern "C" double t_percentile_c(double p, int df);

//==============================================================================
//! Calculate the n-th order Legendre polynomials at the value of x.
//!
//! @param n   The maximum order requested
//! @param x   The value to evaluate at; x is expected to be within [-1,1]
//! @param pnx The requested Legendre polynomials of order 0 to n (inclusive)
//!   evaluated at x.
//==============================================================================

extern "C" void calc_pn_c(int n, double x, double pnx[]);

//==============================================================================
//! Find the value of f(x) given a set of Legendre coefficients and the value
//! of x.
//!
//! @param n    The maximum order of the expansion
//! @param data The polynomial expansion coefficient data; without the (2l+1)/2
//!   factor.
//! @param x    The value to evaluate at; x is expected to be within [-1,1]
//! @return The requested Legendre polynomials of order 0 to n (inclusive)
//!   evaluated at x
//==============================================================================

extern "C" double evaluate_legendre_c(int n, const double data[], double x);

//==============================================================================
//! Calculate the n-th order real spherical harmonics for a given angle (in
//! terms of (u,v,w)) for all 0<=n and -m<=n<=n.
//!
//! @param n      The maximum order requested
//! @param uvw[3] The direction the harmonics are requested at
//! @param rn     The requested harmonics of order 0 to n (inclusive)
//!   evaluated at uvw.
//==============================================================================

extern "C" void calc_rn_c(int n, const double uvw[3], double rn[]);

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
//! @param n   The maximum order requested
//! @param rho The radial parameter to specify location on the unit disk
//! @param phi The angle parameter to specify location on the unit disk
//! @param zn  The requested moments of order 0 to n (inclusive)
//!   evaluated at rho and phi.
//==============================================================================

extern "C" void calc_zn_c(int n, double rho, double phi, double zn[]);

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
//! @param n       The maximum order requested
//! @param rho     The radial parameter to specify location on the unit disk
//! @param phi     The angle parameter to specify location on the unit disk
//! @param zn_rad  The requested moments of order 0 to n (inclusive)
//!   evaluated at rho and phi when m = 0.
//==============================================================================

extern "C" void calc_zn_rad_c(int n, double rho, double zn_rad[]);

//==============================================================================
//! Rotate the direction cosines through a polar angle whose cosine is mu and
//! through an azimuthal angle sampled uniformly.
//!
//! This is done with direct sampling rather than rejection sampling as is done
//! in MCNP and Serpent.
//!
//! @param uvw[3] The initial, and final, direction vector
//! @param mu     The cosine of angle in lab or CM
//! @param phi    The azimuthal angle; will randomly chosen angle if a nullptr
//!   is passed
//==============================================================================

extern "C" void rotate_angle_c(double uvw[3], double mu, double* phi);

Direction rotate_angle(Direction u, double mu, double* phi);

//==============================================================================
//! Samples an energy from the Maxwell fission distribution based on a direct
//! sampling scheme.
//!
//! The probability distribution function for a Maxwellian is given as
//! p(x) = 2/(T*sqrt(pi))*sqrt(x/T)*exp(-x/T). This PDF can be sampled using
//! rule C64 in the Monte Carlo Sampler LA-9721-MS.
//!
//! @param T The tabulated function of the incoming energy
//! @result The sampled outgoing energy
//==============================================================================

extern "C" double maxwell_spectrum_c(double T);

//==============================================================================
//! Samples an energy from a Watt energy-dependent fission distribution.
//!
//! Although fitted parameters exist for many nuclides, generally the
//! continuous tabular distributions (LAW 4) should be used in lieu of the Watt
//! spectrum. This direct sampling scheme is an unpublished scheme based on the
//! original Watt spectrum derivation (See F. Brown's MC lectures).
//!
//! @param a Watt parameter a
//! @param b Watt parameter b
//! @result The sampled outgoing energy
//==============================================================================

extern "C" double watt_spectrum_c(double a, double b);

//==============================================================================
//! Doppler broadens the windowed multipole curvefit.
//!
//! The curvefit is a polynomial of the form a/E + b/sqrt(E) + c + d sqrt(E)...
//!
//! @param E       The energy to evaluate the broadening at
//! @param dopp    sqrt(atomic weight ratio / kT) with kT given in eV
//! @param n       The number of components to the polynomial
//! @param factors The output leading coefficient
//==============================================================================

extern "C" void broaden_wmp_polynomials_c(double E, double dopp, int n,
                                          double factors[]);

//==============================================================================
//! Constructs a natural cubic spline.
//!
//! Given a tabulated function y_i = f(x_i), this computes the second
//! derivative of the interpolating function at each x_i, which can then be
//! used in any subsequent calls to spline_interpolate or spline_integrate for
//! the same set of x and y values.
//!
//! @param n       Number of points
//! @param x       Values of the independent variable, which must be strictly
//!   increasing.
//! @param y       Values of the dependent variable.
//! @param[out] z  The second derivative of the interpolating function at each
//!   value of x.
//==============================================================================

extern "C" void spline_c(int n, const double x[], const double y[], double z[]);

//==============================================================================
//! Determine the cubic spline interpolated y-value for a given x-value.
//!
//! @param n     Number of points
//! @param x     Values of the independent variable, which must be strictly
//!   increasing.
//! @param y     Values of the dependent variable.
//! @param z     The second derivative of the interpolating function at each
//!   value of x.
//! @param xint  Point at which to evaluate the cubic spline polynomial
//! @result      Interpolated value
//==============================================================================

extern "C" double spline_interpolate_c(int n, const double x[], const double y[],
                                       const double z[], double xint);

//==============================================================================
//! Evaluate the definite integral of the interpolating cubic spline between
//! the given endpoints.
//!
//! @param n   Number of points
//! @param x   Values of the independent variable, which must be strictly
//!   increasing.
//! @param y   Values of the dependent variable.
//! @param z   The second derivative of the interpolating function at each
//!   value of x.
//! @param xa  Lower limit of integration
//! @param xb  Upper limit of integration
//! @result    Integral
//==============================================================================

extern "C" double spline_integrate_c(int n, const double x[], const double y[],
                                     const double z[], double xa, double xb);

} // namespace openmc
#endif // OPENMC_MATH_FUNCTIONS_H
