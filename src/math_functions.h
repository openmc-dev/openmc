#ifndef MATH_FUNCTIONS_H
#define MATH_FUNCTIONS_H

#include <cmath>
#include <complex>
#include <iostream>

#include "random_lcg.h"
// #include "faddeeva/Faddeeva.hh"


namespace openmc {


//==============================================================================
// NORMAL_PERCENTILE calculates the percentile of the standard normal
// distribution with a specified probability level
//==============================================================================

extern "C" double normal_percentile_c(double p) __attribute__ ((const));

//==============================================================================
// T_PERCENTILE calculates the percentile of the Student's t distribution with
// a specified probability level and number of degrees of freedom
//==============================================================================

extern "C" double t_percentile_c(double p, int df) __attribute__ ((const));

//==============================================================================
// CALC_PN calculates the n-th order Legendre polynomial at the value of x.
//==============================================================================

extern "C" double calc_pn_c(int n, double x) __attribute__ ((const));

//==============================================================================
// CALC_RN calculates the n-th order spherical harmonics for a given angle
// (in terms of (u,v,w)).  All Rn,m values are provided (where -n<=m<=n)
//==============================================================================

extern "C" void calc_rn_c(int n, double uvw[3], double rn[]);

//==============================================================================
// EVALUATE_LEGENDRE Find the value of f(x) given a set of Legendre coefficients
// and the value of x
//==============================================================================

extern "C" double evaluate_legendre_c(int n, double data[], double x)
     __attribute__ ((const));

//==============================================================================
// ROTATE_ANGLE rotates direction cosines through a polar angle whose cosine is
// mu and through an azimuthal angle sampled uniformly. Note that this is done
// with direct sampling rather than rejection as is done in MCNP and SERPENT.
//==============================================================================

extern "C" void rotate_angle_c(double uvw[3], double mu, double phi = -10.);

//==============================================================================
// MAXWELL_SPECTRUM samples an energy from the Maxwell fission distribution
// based on a direct sampling scheme. The probability distribution function for
// a Maxwellian is given as p(x) = 2/(T*sqrt(pi))*sqrt(x/T)*exp(-x/T).
// This PDF can be sampled using rule C64 in the Monte Carlo Sampler LA-9721-MS.
//==============================================================================

extern "C" double maxwell_spectrum_c(double T);

//==============================================================================
// WATT_SPECTRUM samples the outgoing energy from a Watt energy-dependent
// fission spectrum. Although fitted parameters exist for many nuclides,
// generally the continuous tabular distributions (LAW 4) should be used in
// lieu of the Watt spectrum. This direct sampling scheme is an unpublished
// scheme based on the original Watt spectrum derivation (See F. Brown's
// MC lectures).
//==============================================================================

extern "C" double watt_spectrum_c(double a, double b);

//==============================================================================
// FADDEEVA the Faddeeva function, using Stephen Johnson's implementation
//==============================================================================

// extern "C" std::complex<double> faddeeva_c(std::complex<double> z);

// extern "C" std::complex<double> w_derivative_c(std::complex<double> z, int order);

//==============================================================================
// BROADEN_WMP_POLYNOMIALS Doppler broadens the windowed multipole curvefit.
// The curvefit is a polynomial of the form a/E + b/sqrt(E) + c + d sqrt(E) ...
//==============================================================================

extern "C" void broaden_wmp_polynomials_c(double E, double dopp, int n, double factors[]);

} // namespace openmc
#endif // MATH_FUNCTIONS_H