#ifndef MATH_FUNCTIONS_H
#define MATH_FUNCTIONS_H

#include <cmath>
#include <complex.h>
#include <iostream>

#include "random_lcg.h"
// #include "faddeeva/Faddeeva.h"


namespace openmc {

// TODO: cmath::M_PI has 3 more digits precision than the Fortran constant we
// use so for now we will reuse the Fortran constant until we are OK with
// modifying test results
const double PI = 3.1415926535898;

const double SQRT_PI = std::sqrt(PI);

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
// EVALUATE_LEGENDRE Find the value of f(x) given a set of Legendre coefficients
// and the value of x
//==============================================================================

extern "C" double evaluate_legendre_c(int n, double data[], double x)
     __attribute__ ((const));

//==============================================================================
// CALC_RN calculates the n-th order spherical harmonics for a given angle
// (in terms of (u,v,w)).  All Rn,m values are provided (where -n<=m<=n)
//==============================================================================

extern "C" void calc_rn_c(int n, double uvw[3], double rn[]);

//==============================================================================
// CALC_ZN calculates the n-th order modified Zernike polynomial moment for a
// given angle (rho, theta) location in the unit disk. The normalization of the
// polynomials is such tha the integral of Z_pq*Z_pq over the unit disk is
// exactly pi
//==============================================================================

extern "C" void calc_zn_c(int n, double rho, double phi, double zn[]);

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

// extern "C" double complex faddeeva_c(double complex z) __attribute__ ((const));

// extern "C" double complex w_derivative_c(double complex z, int order);

//==============================================================================
// BROADEN_WMP_POLYNOMIALS Doppler broadens the windowed multipole curvefit.
// The curvefit is a polynomial of the form a/E + b/sqrt(E) + c + d sqrt(E) ...
//==============================================================================

extern "C" void broaden_wmp_polynomials_c(double E, double dopp, int n,
                                          double factors[]);

} // namespace openmc
#endif // MATH_FUNCTIONS_H