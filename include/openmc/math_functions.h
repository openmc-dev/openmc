//! \file math_functions.h
//! A collection of elementary math functions.

#ifndef OPENMC_MATH_FUNCTIONS_H
#define OPENMC_MATH_FUNCTIONS_H

#include <cmath>
#include <complex>
#include <cstdlib>

#include "openmc/constants.h"
#include "openmc/position.h"
#include "openmc/random_lcg.h"


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

extern "C" void rotate_angle_c(double uvw[3], double mu, const double* phi,
  uint64_t* seed);

Direction rotate_angle(Direction u, double mu, const double* phi,
  uint64_t* seed);

//==============================================================================
//! Samples an energy from the Maxwell fission distribution based on a direct
//! sampling scheme.
//!
//! The probability distribution function for a Maxwellian is given as
//! p(x) = 2/(T*sqrt(pi))*sqrt(x/T)*exp(-x/T). This PDF can be sampled using
//! rule C64 in the Monte Carlo Sampler LA-9721-MS.
//!
//! \param T The tabulated function of the incoming energy
//! \param seed A pointer to the pseudorandom seed
//! \return The sampled outgoing energy
//==============================================================================

extern "C" double maxwell_spectrum(double T, uint64_t* seed);

//==============================================================================
//! Samples an energy from a Watt energy-dependent fission distribution.
//!
//! Although fitted parameters exist for many nuclides, generally the
//! continuous tabular distributions (LAW 4) should be used in lieu of the Watt
//! spectrum. This direct sampling scheme is an unpublished scheme based on the
//! original Watt spectrum derivation (See F. Brown's MC lectures).
//!
//! \param a Watt parameter a
//! \param b Watt parameter b
//! \param seed A pointer to the pseudorandom seed
//! \return The sampled outgoing energy
//==============================================================================

extern "C" double watt_spectrum(double a, double b, uint64_t* seed);

//==============================================================================
//! Samples an energy from the Gaussian energy-dependent fission distribution.
//!
//! Samples from a Normal distribution with a given mean and standard deviation
//! The PDF is defined as s(x) = (1/2*sigma*sqrt(2) * e-((mu-x)/2*sigma)^2
//! Its sampled according to
//! http://www-pdg.lbl.gov/2009/reviews/rpp2009-rev-monte-carlo-techniques.pdf
//! section 33.4.4
//!
//! @param mean mean of the Gaussian distribution
//! @param std_dev standard deviation of the Gaussian distribution
//! @param seed A pointer to the pseudorandom seed
//! @result The sampled outgoing energy
//==============================================================================

extern "C" double normal_variate(double mean, double std_dev, uint64_t* seed);

//==============================================================================
//! Samples an energy from the Muir (Gaussian) energy-dependent distribution.
//!
//! This is another form of the Gaussian distribution but with more easily
//! modifiable parameters
//! https://permalink.lanl.gov/object/tr?what=info:lanl-repo/lareport/LA-05411-MS
//!
//! @param e0 peak neutron energy [eV]
//! @param m_rat ratio of the fusion reactants to AMU
//! @param kt the ion temperature of the reactants [eV]
//! @param seed A pointer to the pseudorandom seed
//! @result The sampled outgoing energy
//==============================================================================

extern "C" double muir_spectrum(double e0, double m_rat, double kt,
  uint64_t* seed);

//==============================================================================
//! Doppler broadens the windowed multipole curvefit.
//!
//! The curvefit is a polynomial of the form a/E + b/sqrt(E) + c + d sqrt(E)...
//!
//! \param E       The energy to evaluate the broadening at
//! \param dopp    sqrt(atomic weight ratio / kT) with kT given in eV
//! \param n       The number of components to the polynomial
//! \param factors The output leading coefficient
//==============================================================================

extern "C" void broaden_wmp_polynomials(double E, double dopp, int n, double factors[]);

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

double spline_interpolate(int n, const double x[], const double y[],
  const double z[], double xint);

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

} // namespace openmc
#endif // OPENMC_MATH_FUNCTIONS_H
