#ifndef MATH_FUNCTIONS_H
#define MATH_FUNCTIONS_H

#include <cmath>
#include <cstdlib>

#include "random_lcg.h"


namespace openmc {

//==============================================================================
// Module constants.
//==============================================================================

// TODO: cmath::M_PI has 3 more digits precision than the Fortran constant we
// use so for now we will reuse the Fortran constant until we are OK with
// modifying test results
extern "C" constexpr double PI {3.1415926535898};

extern "C" constexpr double SQRT_PI {std::sqrt(PI)};

//==============================================================================
//! Calculate the percentile of the standard normal distribution with a
//! specified probability level.
//!
//! @param p The probability level
//! @return The requested percentile
//==============================================================================

extern "C" double normal_percentile_c(const double p);

//==============================================================================
//! Calculate the percentile of the Student's t distribution with a specified
//! probability level and number of degrees of freedom.
//!
//! @param p  The probability level
//! @param df The degrees of freedom
//! @return The requested percentile
//==============================================================================

extern "C" double t_percentile_c(const double p, const int df);

//==============================================================================
//! Calculate the n-th order Legendre polynomials at the value of x.
//!
//! @param n   The maximum order requested
//! @param x   The value to evaluate at; x is expected to be within [-1,1]
//! @param pnx The requested Legendre polynomials of order 0 to n (inclusive)
//!   evaluated at x.
//==============================================================================

extern "C" void calc_pn_c(const int n, const double x, double pnx[]);

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

extern "C" double evaluate_legendre_c(const int n, const double data[],
                                      const double x);

//==============================================================================
//! Calculate the n-th order real spherical harmonics for a given angle (in
//! terms of (u,v,w)) for all 0<=n and -m<=n<=n.
//!
//! @param n      The maximum order requested
//! @param uvw[3] The direction the harmonics are requested at
//! @param rn     The requested harmonics of order 0 to n (inclusive)
//!   evaluated at uvw.
//==============================================================================

extern "C" void calc_rn_c(const int n, const double uvw[3], double rn[]);

//==============================================================================
//! Calculate the n-th order modified Zernike polynomial moment for a given
//! angle (rho, theta) location on the unit disk.
//!
//! The normalization of the polynomials is such that the integral of Z_pq^2
//! over the unit disk is exactly pi
//!
//! @param n   The maximum order requested
//! @param rho The radial parameter to specify location on the unit disk
//! @param phi The angle parameter to specify location on the unit disk
//! @param zn  The requested moments of order 0 to n (inclusive)
//!   evaluated at rho and phi.
//==============================================================================

extern "C" void calc_zn_c(const int n, const double rho, const double phi,
                          double zn[]);

//==============================================================================
//! Rotate the direction cosines through a polar angle whose cosine is mu and
//! through an azimuthal angle sampled uniformly.
//!
//! This is done with direct sampling rather than rejection sampling as is done
//! in MCNP and Serpent.
//!
//! @param uvw[3] The initial, and final, direction vector
//! @param mu     The cosine of angle in lab or CM
//! @param phi    The azimuthal angle; defaults to a randomly chosen angle
//==============================================================================

extern "C" void rotate_angle_c(double uvw[3], const double mu,
                               double* phi=nullptr);

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

extern "C" double maxwell_spectrum_c(const double T);

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

extern "C" double watt_spectrum_c(const double a, const double b);

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

extern "C" void broaden_wmp_polynomials_c(const double E, const double dopp,
                                          const int n, double factors[]);

} // namespace openmc
#endif // MATH_FUNCTIONS_H