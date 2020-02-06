#include "openmc/math_functions.h"

#include "Faddeeva.hh"

namespace openmc {

//==============================================================================
// Mathematical methods
//==============================================================================

double normal_percentile(double p)
{
  constexpr double p_low = 0.02425;
  constexpr double a[6] = {-3.969683028665376e1, 2.209460984245205e2,
                           -2.759285104469687e2, 1.383577518672690e2,
                           -3.066479806614716e1, 2.506628277459239e0};
  constexpr double b[5] = {-5.447609879822406e1, 1.615858368580409e2,
                           -1.556989798598866e2, 6.680131188771972e1,
                           -1.328068155288572e1};
  constexpr double c[6] = {-7.784894002430293e-3, -3.223964580411365e-1,
                           -2.400758277161838, -2.549732539343734,
                           4.374664141464968, 2.938163982698783};
  constexpr double d[4] = {7.784695709041462e-3, 3.224671290700398e-1,
                           2.445134137142996, 3.754408661907416};

  // The rational approximation used here is from an unpublished work at
  // http://home.online.no/~pjacklam/notes/invnorm/

  double z;
  double q;

  if (p < p_low) {
    // Rational approximation for lower region.

    q = std::sqrt(-2.0 * std::log(p));
    z = (((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) /
          ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1.0);

  } else if (p <= 1.0 - p_low) {
    // Rational approximation for central region
    q = p - 0.5;
    double r = q * q;
    z = (((((a[0]*r + a[1])*r + a[2])*r + a[3])*r + a[4])*r + a[5])*q /
         (((((b[0]*r + b[1])*r + b[2])*r + b[3])*r + b[4])*r + 1.0);

  } else {
    // Rational approximation for upper region

    q = std::sqrt(-2.0*std::log(1.0 - p));
    z = -(((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) /
          ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1.0);
  }

  // Refinement based on Newton's method

  z = z - (0.5 * std::erfc(-z / std::sqrt(2.0)) - p) * std::sqrt(2.0 * PI) *
       std::exp(0.5 * z * z);

  return z;

}


double t_percentile(double p, int df)
{
  double t;

  if (df == 1) {
    // For one degree of freedom, the t-distribution becomes a Cauchy
    // distribution whose cdf we can invert directly

    t = std::tan(PI*(p - 0.5));
  } else if (df == 2) {
     // For two degrees of freedom, the cdf is given by 1/2 + x/(2*sqrt(x^2 +
     // 2)). This can be directly inverted to yield the solution below

    t = 2.0 * std::sqrt(2.0)*(p - 0.5) /
         std::sqrt(1. - 4. * std::pow(p - 0.5, 2.));
  } else {
    // This approximation is from E. Olusegun George and Meenakshi Sivaram, "A
    // modification of the Fisher-Cornish approximation for the student t
    // percentiles," Communication in Statistics - Simulation and Computation,
    // 16 (4), pp. 1123-1132 (1987).
    double n = df;
    double k = 1. / (n - 2.);
    double z = normal_percentile(p);
    double z2 = z * z;
    t = std::sqrt(n * k) * (z + (z2 - 3.) * z * k / 4. + ((5. * z2 - 56.) * z2 +
         75.) * z * k * k / 96. + (((z2 - 27.) * 3. * z2 + 417.) * z2 - 315.) *
         z * k * k * k / 384.);
  }

  return t;
}


void calc_pn_c(int n, double x, double pnx[])
{
  pnx[0] = 1.;
  if (n >= 1) {
    pnx[1] = x;
  }

  // Use recursion relation to build the higher orders
  for (int l = 1; l < n; l++) {
    pnx[l + 1] = ((2 * l + 1) * x * pnx[l] - l * pnx[l - 1]) / (l + 1);
  }
}


double evaluate_legendre(int n, const double data[], double x)
{
  double* pnx = new double[n + 1];
  double val = 0.0;
  calc_pn_c(n, x, pnx);
  for (int l = 0; l <= n; l++) {
    val += (l + 0.5) * data[l] * pnx[l];
  }
  delete[] pnx;
  return val;
}


void calc_rn_c(int n, const double uvw[3], double rn[])
{
  Direction u {uvw};
  calc_rn(n, u, rn);
}


void calc_rn(int n, Direction u, double rn[])
{
  // rn[] is assumed to have already been allocated to the correct size

  // Store the cosine of the polar angle and the azimuthal angle
  double w = u.z;
  double phi;
  if (u.x == 0.) {
    phi = 0.;
  } else {
    phi = std::atan2(u.y, u.x);
  }

  // Store the shorthand of 1-w * w
  double w2m1 = 1. - w * w;

  // Now evaluate the spherical harmonics function
  rn[0] = 1.;
  int i = 0;
  for (int l = 1; l <= n; l++) {
    // Set the index to the start of this order
    i += 2 * (l - 1) + 1;

    // Now evaluate each
    switch (l) {
      case 1:
        // l = 1, m = -1
        rn[i] = -(std::sqrt(w2m1) * std::sin(phi));
        // l = 1, m = 0
        rn[i + 1] = w;
        // l = 1, m = 1
        rn[i + 2] = -(std::sqrt(w2m1) * std::cos(phi));
        break;
      case 2:
        // l = 2, m = -2
        rn[i] = 0.288675134594813 * (-3. * w * w + 3.) * std::sin(2. * phi);
        // l = 2, m = -1
        rn[i + 1] = -(1.73205080756888 * w*std::sqrt(w2m1) * std::sin(phi));
        // l = 2, m = 0
        rn[i + 2] = 1.5 * w * w - 0.5;
        // l = 2, m = 1
        rn[i + 3] = -(1.73205080756888 * w*std::sqrt(w2m1) * std::cos(phi));
        // l = 2, m = 2
        rn[i + 4] = 0.288675134594813 * (-3. * w * w + 3.) * std::cos(2. * phi);
        break;
      case 3:
        // l = 3, m = -3
        rn[i] = -(0.790569415042095 * std::pow(w2m1, 1.5) * std::sin(3. * phi));
        // l = 3, m = -2
        rn[i + 1] = 1.93649167310371 * w*(w2m1) * std::sin(2.*phi);
        // l = 3, m = -1
        rn[i + 2] = -(0.408248290463863*std::sqrt(w2m1)*((7.5)*w * w - 3./2.) *
             std::sin(phi));
        // l = 3, m = 0
        rn[i + 3] = 2.5 * std::pow(w, 3) - 1.5 * w;
        // l = 3, m = 1
        rn[i + 4] = -(0.408248290463863*std::sqrt(w2m1)*((7.5)*w * w - 3./2.) *
             std::cos(phi));
        // l = 3, m = 2
        rn[i + 5] = 1.93649167310371 * w*(w2m1) * std::cos(2.*phi);
        // l = 3, m = 3
        rn[i + 6] = -(0.790569415042095 * std::pow(w2m1, 1.5) * std::cos(3.* phi));
        break;
      case 4:
        // l = 4, m = -4
        rn[i] = 0.739509972887452 * (w2m1 * w2m1) * std::sin(4.0*phi);
        // l = 4, m = -3
        rn[i + 1] = -(2.09165006633519 * w * std::pow(w2m1, 1.5) * std::sin(3.* phi));
        // l = 4, m = -2
        rn[i + 2] = 0.074535599249993 * (w2m1)*(52.5 * w * w - 7.5) * std::sin(2. *phi);
        // l = 4, m = -1
        rn[i + 3] = -(0.316227766016838*std::sqrt(w2m1)*(17.5 * std::pow(w, 3) - 7.5 * w) *
             std::sin(phi));
        // l = 4, m = 0
        rn[i + 4] = 4.375 * std::pow(w, 4) - 3.75 * w * w + 0.375;
        // l = 4, m = 1
        rn[i + 5] = -(0.316227766016838*std::sqrt(w2m1)*(17.5 * std::pow(w, 3) - 7.5*w) *
             std::cos(phi));
        // l = 4, m = 2
        rn[i + 6] = 0.074535599249993 * (w2m1)*(52.5*w * w - 7.5) * std::cos(2.*phi);
        // l = 4, m = 3
        rn[i + 7] = -(2.09165006633519 * w * std::pow(w2m1, 1.5) * std::cos(3.* phi));
        // l = 4, m = 4
        rn[i + 8] = 0.739509972887452 * w2m1 * w2m1 * std::cos(4.0*phi);
        break;
      case 5:
        // l = 5, m = -5
        rn[i] = -(0.701560760020114 * std::pow(w2m1, 2.5) * std::sin(5.0 * phi));
        // l = 5, m = -4
        rn[i + 1] = 2.21852991866236 * w * w2m1 * w2m1 * std::sin(4.0 * phi);
        // l = 5, m = -3
        rn[i + 2] = -(0.00996023841111995 * std::pow(w2m1, 1.5) *
             ((945.0 /2.)* w * w - 52.5) * std::sin(3.*phi));
        // l = 5, m = -2
        rn[i + 3] = 0.0487950036474267 * (w2m1)
             * ((315.0/2.)* std::pow(w, 3) - 52.5 * w) * std::sin(2.*phi);
        // l = 5, m = -1
        rn[i + 4] = -(0.258198889747161*std::sqrt(w2m1) *
             (39.375 * std::pow(w, 4) - 105.0/4.0 * w * w + 15.0/8.0) * std::sin(phi));
        // l = 5, m = 0
        rn[i + 5] = 7.875 * std::pow(w, 5) - 8.75 * std::pow(w, 3) + 1.875 * w;
        // l = 5, m = 1
        rn[i + 6] = -(0.258198889747161 * std::sqrt(w2m1)*
             (39.375 * std::pow(w, 4) - 105.0/4.0 * w * w + 15.0/8.0) * std::cos(phi));
        // l = 5, m = 2
        rn[i + 7] = 0.0487950036474267 * (w2m1) *
             ((315.0 / 2.) * std::pow(w, 3) - 52.5*w) * std::cos(2.*phi);
        // l = 5, m = 3
        rn[i + 8] = -(0.00996023841111995 * std::pow(w2m1, 1.5) *
             ((945.0 / 2.) * w * w - 52.5) * std::cos(3.*phi));
        // l = 5, m = 4
        rn[i + 9] = 2.21852991866236 * w * w2m1 * w2m1 * std::cos(4.0*phi);
        // l = 5, m = 5
        rn[i + 10] = -(0.701560760020114 * std::pow(w2m1, 2.5) * std::cos(5.0* phi));
        break;
      case 6:
        // l = 6, m = -6
        rn[i] = 0.671693289381396 * std::pow(w2m1, 3) * std::sin(6.0*phi);
        // l = 6, m = -5
        rn[i + 1] = -(2.32681380862329 * w*std::pow(w2m1, 2.5) * std::sin(5.0*phi));
        // l = 6, m = -4
        rn[i + 2] = 0.00104990131391452 * w2m1 * w2m1 *
             ((10395.0/2.) * w * w - 945.0/2.) * std::sin(4.0 * phi);
        // l = 6, m = -3
        rn[i + 3] = -(0.00575054632785295 * std::pow(w2m1, 1.5) *
             ((3465.0/2.) * std::pow(w, 3) - 945.0/2.*w) * std::sin(3.*phi));
        // l = 6, m = -2
        rn[i + 4] = 0.0345032779671177 * (w2m1) *
             ((3465.0/8.0)* std::pow(w, 4) - 945.0/4.0 * w * w + 105.0/8.0) *
             std::sin(2. * phi);
        // l = 6, m = -1
        rn[i + 5] = -(0.218217890235992*std::sqrt(w2m1) *
             ((693.0/8.0)* std::pow(w, 5)- 315.0/4.0 * std::pow(w, 3) + (105.0/8.0)*w) *
             std::sin(phi));
        // l = 6, m = 0
        rn[i + 6] = 14.4375 * std::pow(w, 6) - 19.6875 * std::pow(w, 4) + 6.5625 * w * w -
             0.3125;
        // l = 6, m = 1
        rn[i + 7] = -(0.218217890235992*std::sqrt(w2m1) *
             ((693.0/8.0)* std::pow(w, 5)- 315.0/4.0 * std::pow(w, 3) + (105.0/8.0)*w) *
             std::cos(phi));
        // l = 6, m = 2
        rn[i + 8] = 0.0345032779671177 * w2m1 *
             ((3465.0/8.0)* std::pow(w, 4) -945.0/4.0 * w * w + 105.0/8.0) *
             std::cos(2.*phi);
        // l = 6, m = 3
        rn[i + 9] = -(0.00575054632785295 * std::pow(w2m1, 1.5) *
             ((3465.0/2.) * std::pow(w, 3) - 945.0/2.*w) * std::cos(3.*phi));
        // l = 6, m = 4
        rn[i + 10] = 0.00104990131391452 * w2m1 * w2m1 *
             ((10395.0/2.)*w * w - 945.0/2.) * std::cos(4.0*phi);
        // l = 6, m = 5
        rn[i + 11] = -(2.32681380862329 * w * std::pow(w2m1, 2.5) * std::cos(5.0*phi));
        // l = 6, m = 6
        rn[i + 12] = 0.671693289381396 * std::pow(w2m1, 3) * std::cos(6.0*phi);
        break;
      case 7:
        // l = 7, m = -7
        rn[i] = -(0.647259849287749 * std::pow(w2m1, 3.5) * std::sin(7.0*phi));
        // l = 7, m = -6
        rn[i + 1] = 2.42182459624969 * w*std::pow(w2m1, 3) * std::sin(6.0*phi);
        // l = 7, m = -5
        rn[i + 2] = -(9.13821798555235e-5*std::pow(w2m1, 2.5) *
             ((135135.0/2.)*w * w - 10395.0/2.) * std::sin(5.0*phi));
        // l = 7, m = -4
        rn[i + 3] = 0.000548293079133141 * w2m1 * w2m1 *
             ((45045.0/2.)*std::pow(w, 3) - 10395.0/2.*w) * std::sin(4.0*phi);
        // l = 7, m = -3
        rn[i + 4] = -(0.00363696483726654 * std::pow(w2m1, 1.5) *
             ((45045.0/8.0)* std::pow(w, 4) - 10395.0/4.0 * w * w + 945.0/8.0) *
             std::sin(3.*phi));
        // l = 7, m = -2
        rn[i + 5] = 0.025717224993682 * (w2m1) *
             ((9009.0/8.0)* std::pow(w, 5) -3465.0/4.0 * std::pow(w, 3) + (945.0/8.0)*w) *
             std::sin(2.*phi);
        // l = 7, m = -1
        rn[i + 6] = -(0.188982236504614*std::sqrt(w2m1) *
             ((3003.0/16.0)* std::pow(w, 6) - 3465.0/16.0 * std::pow(w, 4) +
             (945.0/16.0)*w * w - 35.0/16.0) * std::sin(phi));
        // l = 7, m = 0
        rn[i + 7] = 26.8125 * std::pow(w, 7) - 43.3125 * std::pow(w, 5) + 19.6875 * std::pow(w, 3) -
             2.1875 * w;
        // l = 7, m = 1
        rn[i + 8] = -(0.188982236504614*std::sqrt(w2m1) * ((3003.0/16.0) * std::pow(w, 6) -
             3465.0/16.0 * std::pow(w, 4) + (945.0/16.0)*w * w - 35.0/16.0) * std::cos(phi));
        // l = 7, m = 2
        rn[i + 9] = 0.025717224993682 * (w2m1) * ((9009.0/8.0)* std::pow(w, 5) -
             3465.0/4.0 * std::pow(w, 3) + (945.0/8.0)*w) * std::cos(2.*phi);
        // l = 7, m = 3
        rn[i + 10] = -(0.00363696483726654 * std::pow(w2m1, 1.5) *
             ((45045.0/8.0)* std::pow(w, 4) - 10395.0/4.0 * w * w + 945.0/8.0) *
             std::cos(3.*phi));
        // l = 7, m = 4
        rn[i + 11] = 0.000548293079133141 * w2m1 * w2m1 *
             ((45045.0/2.)*std::pow(w, 3) - 10395.0/2.*w) * std::cos(4.0*phi);
        // l = 7, m = 5
        rn[i + 12] = -(9.13821798555235e-5*std::pow(w2m1, 2.5) *
             ((135135.0/2.)*w * w - 10395.0/2.) * std::cos(5.0*phi));
        // l = 7, m = 6
        rn[i + 13] = 2.42182459624969 * w*std::pow(w2m1, 3) * std::cos(6.0*phi);
        // l = 7, m = 7
        rn[i + 14] = -(0.647259849287749 * std::pow(w2m1, 3.5) * std::cos(7.0*phi));
        break;
      case 8:
        // l = 8, m = -8
        rn[i] = 0.626706654240044 * std::pow(w2m1, 4) * std::sin(8.0*phi);
        // l = 8, m = -7
        rn[i + 1] = -(2.50682661696018 * w*std::pow(w2m1, 3.5) * std::sin(7.0*phi));
        // l = 8, m = -6
        rn[i + 2] = 6.77369783729086e-6*std::pow(w2m1, 3)*
             ((2027025.0/2.)*w * w - 135135.0/2.) * std::sin(6.0*phi);
        // l = 8, m = -5
        rn[i + 3] = -(4.38985792528482e-5*std::pow(w2m1, 2.5) *
                  ((675675.0/2.)*std::pow(w, 3) - 135135.0/2.*w) * std::sin(5.0*phi));
        // l = 8, m = -4
        rn[i + 4] = 0.000316557156832328 * w2m1 * w2m1 *
             ((675675.0/8.0)* std::pow(w, 4) - 135135.0/4.0 * w * w + 10395.0/8.0) *
             std::sin(4.0*phi);
        // l = 8, m = -3
        rn[i + 5] = -(0.00245204119306875 * std::pow(w2m1, 1.5) * ((135135.0/8.0) *
             std::pow(w, 5) - 45045.0/4.0 * std::pow(w, 3) + (10395.0/8.0)*w) * std::sin(3.*phi));
        // l = 8, m = -2
        rn[i + 6] = 0.0199204768222399 * (w2m1) *
             ((45045.0/16.0)* std::pow(w, 6)- 45045.0/16.0 * std::pow(w, 4) +
             (10395.0/16.0)*w * w - 315.0/16.0) * std::sin(2.*phi);
        // l = 8, m = -1
        rn[i + 7] = -(0.166666666666667*std::sqrt(w2m1) *
             ((6435.0/16.0)* std::pow(w, 7) - 9009.0/16.0 * std::pow(w, 5) +
             (3465.0/16.0)*std::pow(w, 3) - 315.0/16.0 * w) * std::sin(phi));
        // l = 8, m = 0
        rn[i + 8] = 50.2734375 * std::pow(w, 8) - 93.84375 * std::pow(w, 6) + 54.140625 *
             std::pow(w, 4) - 9.84375 * w * w + 0.2734375;
        // l = 8, m = 1
        rn[i + 9] = -(0.166666666666667*std::sqrt(w2m1) *
                   ((6435.0/16.0)* std::pow(w, 7) - 9009.0/16.0 * std::pow(w, 5) +
                   (3465.0/16.0)*std::pow(w, 3) - 315.0/16.0 * w) * std::cos(phi));
        // l = 8, m = 2
        rn[i + 10] = 0.0199204768222399 * (w2m1)*((45045.0/16.0)* std::pow(w, 6)-
             45045.0/16.0 * std::pow(w, 4) + (10395.0/16.0)*w * w -
             315.0/16.0) * std::cos(2.*phi);
        // l = 8, m = 3
        rn[i + 11] = -(0.00245204119306875 * std::pow(w2m1, 1.5)*
                   ((135135.0/8.0) * std::pow(w, 5) - 45045.0/4.0 * std::pow(w, 3) +
                   (10395.0/8.0)*w) * std::cos(3.*phi));
        // l = 8, m = 4
        rn[i + 12] = 0.000316557156832328 * w2m1 * w2m1*((675675.0/8.0)* std::pow(w, 4) -
             135135.0/4.0 * w * w + 10395.0/8.0) * std::cos(4.0*phi);
        // l = 8, m = 5
        rn[i + 13] = -(4.38985792528482e-5*std::pow(w2m1, 2.5)*((675675.0/2.)*std::pow(w, 3) -
                   135135.0/2.*w) * std::cos(5.0*phi));
        // l = 8, m = 6
        rn[i + 14] = 6.77369783729086e-6*std::pow(w2m1, 3)*((2027025.0/2.)*w * w -
             135135.0/2.) * std::cos(6.0*phi);
        // l = 8, m = 7
        rn[i + 15] = -(2.50682661696018 * w*std::pow(w2m1, 3.5) * std::cos(7.0*phi));
        // l = 8, m = 8
        rn[i + 16] = 0.626706654240044 * std::pow(w2m1, 4) * std::cos(8.0*phi);
        break;
      case 9:
        // l = 9, m = -9
        rn[i] = -(0.609049392175524 * std::pow(w2m1, 4.5) * std::sin(9.0 * phi));
        // l = 9, m = -8
        rn[i + 1] = 2.58397773170915 * w*std::pow(w2m1, 4) * std::sin(8.0 * phi);
        // l = 9, m = -7
        rn[i + 2] = -(4.37240315267812e-7*std::pow(w2m1, 3.5) *
             ((34459425.0/2.)*w * w - 2027025.0/2.) * std::sin(7.0 * phi));
        // l = 9, m = -6
        rn[i + 3] = 3.02928976464514e-6*std::pow(w2m1, 3)*
             ((11486475.0/2.)*std::pow(w, 3) - 2027025.0/2.*w) * std::sin(6.0 * phi);
        // l = 9, m = -5
        rn[i + 4] = -(2.34647776186144e-5*std::pow(w2m1, 2.5) *
             ((11486475.0/8.0)* std::pow(w, 4) - 2027025.0 / 4.0 * w * w +
              135135.0/8.0) * std::sin(5.0 * phi));
        // l = 9, m = -4
        rn[i + 5] = 0.000196320414650061 * w2m1 * w2m1*((2297295.0/8.0)* std::pow(w, 5) -
             675675.0/4.0 * std::pow(w, 3) + (135135.0/8.0)*w) * std::sin(4.0*phi);
        // l = 9, m = -3
        rn[i + 6] = -(0.00173385495536766 * std::pow(w2m1, 1.5) *
                  ((765765.0/16.0)* std::pow(w, 6) - 675675.0/16.0 * std::pow(w, 4) +
                  (135135.0/16.0)*w * w - 3465.0/16.0) * std::sin(3. * phi));
        // l = 9, m = -2
        rn[i + 7] = 0.0158910431540932 * (w2m1)*((109395.0/16.0)* std::pow(w, 7)-
             135135.0/16.0 * std::pow(w, 5) + (45045.0/16.0)*std::pow(w, 3) -
             3465.0/16.0 * w) * std::sin(2. * phi);
        // l = 9, m = -1
        rn[i + 8] = -(0.149071198499986*std::sqrt(w2m1)*((109395.0/128.0)* std::pow(w, 8) -
                  45045.0/32.0 * std::pow(w, 6) + (45045.0/64.0)* std::pow(w, 4) -
                  3465.0/32.0 * w * w + 315.0/128.0) * std::sin(phi));
        // l = 9, m = 0
        rn[i + 9] = 94.9609375 * std::pow(w, 9) - 201.09375 * std::pow(w, 7) +
             140.765625 * std::pow(w, 5)- 36.09375 * std::pow(w, 3) + 2.4609375 * w;
        // l = 9, m = 1
        rn[i + 10] = -(0.149071198499986*std::sqrt(w2m1)*((109395.0/128.0)* std::pow(w, 8) -
                   45045.0/32.0 * std::pow(w, 6) + (45045.0/64.0)* std::pow(w, 4) -
                   3465.0/32.0 * w * w + 315.0/128.0) * std::cos(phi));
        // l = 9, m = 2
        rn[i + 11] = 0.0158910431540932 * (w2m1)*((109395.0/16.0)* std::pow(w, 7) -
             135135.0/16.0 * std::pow(w, 5) + (45045.0/16.0)*std::pow(w, 3) -
             3465.0/ 16.0 * w) * std::cos(2. * phi);
        // l = 9, m = 3
        rn[i + 12] = -(0.00173385495536766 * std::pow(w2m1, 1.5)*((765765.0/16.0) *
                   std::pow(w, 6) - 675675.0/16.0 * std::pow(w, 4) +
                   (135135.0/16.0)* w * w - 3465.0/16.0)* std::cos(3. * phi));
        // l = 9, m = 4
        rn[i + 13] = 0.000196320414650061 * w2m1 * w2m1*((2297295.0/8.0) * std::pow(w, 5) -
             675675.0/4.0 * std::pow(w, 3) + (135135.0/8.0)*w) * std::cos(4.0 * phi);
        // l = 9, m = 5
        rn[i + 14] = -(2.34647776186144e-5*std::pow(w2m1, 2.5)*((11486475.0/8.0) *
                   std::pow(w, 4) - 2027025.0/4.0 * w * w + 135135.0/8.0) *
                   std::cos(5.0 * phi));
        // l = 9, m = 6
        rn[i + 15] = 3.02928976464514e-6*std::pow(w2m1, 3)*((11486475.0/2.)*std::pow(w, 3) -
             2027025.0/2. * w) * std::cos(6.0 * phi);
        // l = 9, m = 7
        rn[i + 16] = -(4.37240315267812e-7*std::pow(w2m1, 3.5)*
                   ((34459425.0/2.) * w * w - 2027025.0/2.) * std::cos(7.0 * phi));
        // l = 9, m = 8
        rn[i + 17] = 2.58397773170915 * w*std::pow(w2m1, 4) * std::cos(8.0 * phi);
        // l = 9, m = 9
        rn[i + 18] = -(0.609049392175524 * std::pow(w2m1, 4.5) * std::cos(9.0 * phi));
        break;
      case 10:
        // l = 10, m = -10
        rn[i] = 0.593627917136573 * std::pow(w2m1, 5) * std::sin(10.0 * phi);
        // l = 10, m = -9
        rn[i + 1] = -(2.65478475211798 * w * std::pow(w2m1, 4.5) * std::sin(9.0 * phi));
        // l = 10, m = -8
        rn[i + 2] = 2.49953651452314e-8 * std::pow(w2m1, 4) *
             ((654729075.0/2.) * w * w - 34459425.0/2.) * std::sin(8.0 * phi);
        // l = 10, m = -7
        rn[i + 3] = -(1.83677671621093e-7*std::pow(w2m1, 3.5)*
                  ((218243025.0/2.)*std::pow(w, 3) - 34459425.0/2.*w) *
                  std::sin(7.0 * phi));
        // l = 10, m = -6
        rn[i + 4] = 1.51464488232257e-6*std::pow(w2m1, 3)*((218243025.0/8.0)* std::pow(w, 4) -
             34459425.0/4.0 * w * w + 2027025.0/8.0) * std::sin(6.0 * phi);
        // l = 10, m = -5
        rn[i + 5] = -(1.35473956745817e-5*std::pow(w2m1, 2.5)*
                  ((43648605.0/8.0)* std::pow(w, 5) - 11486475.0/4.0 * std::pow(w, 3) +
                  (2027025.0/8.0)*w) * std::sin(5.0 * phi));
        // l = 10, m = -4
        rn[i + 6] = 0.000128521880085575 * w2m1 * w2m1*((14549535.0/16.0)* std::pow(w, 6) -
             11486475.0/16.0 * std::pow(w, 4) + (2027025.0/16.0)*w * w -
             45045.0/16.0) * std::sin(4.0 * phi);
        // l = 10, m = -3
        rn[i + 7] = -(0.00127230170115096 * std::pow(w2m1, 1.5)*
                  ((2078505.0/16.0)* std::pow(w, 7) - 2297295.0/16.0 * std::pow(w, 5) +
                  (675675.0/16.0)*std::pow(w, 3) - 45045.0/16.0 * w) * std::sin(3. * phi));
        // l = 10, m = -2
        rn[i + 8] = 0.012974982402692 * (w2m1)*((2078505.0/128.0)* std::pow(w, 8) -
             765765.0/32.0 * std::pow(w, 6) + (675675.0/64.0)* std::pow(w, 4) -
             45045.0/32.0 * w * w + 3465.0/128.0) * std::sin(2. * phi);
        // l = 10, m = -1
        rn[i + 9] = -(0.134839972492648*std::sqrt(w2m1)*((230945.0/128.0)* std::pow(w, 9) -
                   109395.0/32.0 * std::pow(w, 7) + (135135.0/64.0)* std::pow(w, 5) -
                   15015.0/32.0 * std::pow(w, 3) + (3465.0/128.0)*w) * std::sin(phi));
        // l = 10, m = 0
        rn[i + 10] = 180.42578125 * std::pow(w, 10) - 427.32421875 * std::pow(w, 8) +351.9140625
             * std::pow(w, 6) - 117.3046875 * std::pow(w, 4) + 13.53515625 * w * w -0.24609375;
        // l = 10, m = 1
        rn[i + 11] = -(0.134839972492648*std::sqrt(w2m1)*((230945.0/128.0)* std::pow(w, 9) -
                   109395.0/32.0 * std::pow(w, 7) + (135135.0/64.0)* std::pow(w, 5) -15015.0/
                   32.0 * std::pow(w, 3) + (3465.0/128.0)*w) * std::cos(phi));
        // l = 10, m = 2
        rn[i + 12] = 0.012974982402692 * (w2m1)*((2078505.0/128.0)* std::pow(w, 8) -
             765765.0/32.0 * std::pow(w, 6) + (675675.0/64.0)* std::pow(w, 4) -
             45045.0/32.0 * w * w + 3465.0/128.0) * std::cos(2. * phi);
        // l = 10, m = 3
        rn[i + 13] = -(0.00127230170115096 * std::pow(w2m1, 1.5)*
                   ((2078505.0/16.0)* std::pow(w, 7) - 2297295.0/16.0 * std::pow(w, 5) +
                   (675675.0/16.0)*std::pow(w, 3) - 45045.0/16.0 * w) * std::cos(3. * phi));
        // l = 10, m = 4
        rn[i + 14] = 0.000128521880085575 * w2m1 * w2m1*((14549535.0/16.0)* std::pow(w, 6) -
             11486475.0/16.0 * std::pow(w, 4) + (2027025.0/16.0) * w * w -
             45045.0/16.0) * std::cos(4.0 * phi);
        // l = 10, m = 5
        rn[i + 15] = -(1.35473956745817e-5*std::pow(w2m1, 2.5)*
                   ((43648605.0/8.0)* std::pow(w, 5) - 11486475.0/4.0 * std::pow(w, 3) +
                   (2027025.0/8.0)*w) * std::cos(5.0 * phi));
        // l = 10, m = 6
        rn[i + 16] = 1.51464488232257e-6*std::pow(w2m1, 3)*((218243025.0/8.0)* std::pow(w, 4) -
             34459425.0/4.0 * w * w + 2027025.0/8.0) * std::cos(6.0 * phi);
        // l = 10, m = 7
        rn[i + 17] = -(1.83677671621093e-7*std::pow(w2m1, 3.5) *
             ((218243025.0/2.)*std::pow(w, 3) - 34459425.0/2.*w) * std::cos(7.0 * phi));
        // l = 10, m = 8
        rn[i + 18] = 2.49953651452314e-8*std::pow(w2m1, 4)*
             ((654729075.0/2.)*w * w - 34459425.0/2.) * std::cos(8.0 * phi);
        // l = 10, m = 9
        rn[i + 19] = -(2.65478475211798 * w*std::pow(w2m1, 4.5) * std::cos(9.0 * phi));
        // l = 10, m = 10
        rn[i + 20] = 0.593627917136573 * std::pow(w2m1, 5) * std::cos(10.0 * phi);
    }
  }
}


void calc_zn(int n, double rho, double phi, double zn[]) {
  // ===========================================================================
  // Determine vector of sin(n*phi) and cos(n*phi). This takes advantage of the
  // following recurrence relations so that only a single sin/cos have to be
  // evaluated (http://mathworld.wolfram.com/Multiple-AngleFormulas.html)
  //
  // sin(nx) = 2 cos(x) sin((n-1)x) - sin((n-2)x)
  // cos(nx) = 2 cos(x) cos((n-1)x) - cos((n-2)x)

  double sin_phi = std::sin(phi);
  double cos_phi = std::cos(phi);

  std::vector<double> sin_phi_vec(n + 1); // Sin[n * phi]
  std::vector<double> cos_phi_vec(n + 1); // Cos[n * phi]
  sin_phi_vec[0] = 1.0;
  cos_phi_vec[0] = 1.0;
  sin_phi_vec[1] = 2.0 * cos_phi;
  cos_phi_vec[1] = cos_phi;

  for (int i = 2; i <= n; i++) {
    sin_phi_vec[i] = 2. * cos_phi * sin_phi_vec[i - 1] - sin_phi_vec[i - 2];
    cos_phi_vec[i] = 2. * cos_phi * cos_phi_vec[i - 1] - cos_phi_vec[i - 2];
  }

  for (int i = 0; i <= n; i++) {
    sin_phi_vec[i] *=  sin_phi;
  }

  // ===========================================================================
  // Calculate R_pq(rho)
  // Matrix forms of the coefficients which are easier to work with
  std::vector<std::vector<double>> zn_mat(n + 1, std::vector<double>(n + 1));

  // Fill the main diagonal first (Eq 3.9 in Chong)
  for (int p = 0; p <= n; p++) {
    zn_mat[p][p] = std::pow(rho, p);
  }

  // Fill the 2nd diagonal (Eq 3.10 in Chong)
  for (int q = 0; q <= n - 2; q++) {
    zn_mat[q][q+2] = (q + 2) * zn_mat[q+2][q+2] - (q + 1) * zn_mat[q][q];
  }

  // Fill in the rest of the values using the original results (Eq. 3.8 in Chong)
  for (int p = 4; p <= n; p++) {
    double k2 = 2 * p * (p - 1) * (p - 2);
    for (int q = p - 4; q >= 0; q -= 2) {
      double k1 = ((p + q) * (p - q) * (p - 2)) / 2.;
      double k3 = -q * q * (p - 1) - p * (p - 1) * (p - 2);
      double k4 = (-p * (p + q - 2) * (p - q - 2)) / 2.;
      zn_mat[q][p] =
        ((k2 * rho * rho + k3) * zn_mat[q][p-2] + k4 * zn_mat[q][p-4]) / k1;
    }
  }

  // Roll into a single vector for easier computation later
  // The vector is ordered (0,0), (1,-1), (1,1), (2,-2), (2,0),
  // (2, 2), ....   in (n,m) indices
  // Note that the cos and sin vectors are offset by one
  // sin_phi_vec = [sin(x), sin(2x), sin(3x) ...]
  // cos_phi_vec = [1.0, cos(x), cos(2x)... ]
  int i = 0;
  for (int p = 0; p <= n; p++) {
    for (int q = -p; q <= p; q += 2) {
      if (q < 0) {
        zn[i] = zn_mat[std::abs(q)][p] * sin_phi_vec[std::abs(q) - 1];
      } else if (q == 0) {
        zn[i] = zn_mat[q][p];
      } else {
        zn[i] = zn_mat[q][p] * cos_phi_vec[q];
      }
      i++;
    }
  }

}

void calc_zn_rad(int n, double rho, double zn_rad[]) {
  // Calculate R_p0(rho) as Zn_p0(rho)
  // Set up the array of the coefficients

  double q = 0;

  // R_00 is always 1
  zn_rad[0] = 1;

  // Fill in the rest of the array (Eq 3.8 and Eq 3.10 in Chong)
  for (int p = 2; p <= n; p += 2) {
    int index = int(p/2);
    if (p == 2) {
    // Setting up R_22 to calculate R_20 (Eq 3.10 in Chong)
      double R_22 = rho * rho;
      zn_rad[index] = 2 * R_22 - zn_rad[0];
    } else {
      double k1 = ((p + q) * (p - q) * (p - 2)) / 2.;
      double k2 = 2 * p * (p - 1) * (p - 2);
      double k3 = -q * q * (p - 1) - p * (p - 1) * (p - 2);
      double k4 = (-p * (p + q - 2) * (p - q - 2)) / 2.;
      zn_rad[index] =
        ((k2 * rho * rho + k3) * zn_rad[index-1] + k4 * zn_rad[index-2]) / k1;
    }
  }
}


void rotate_angle_c(double uvw[3], double mu, const double* phi, uint64_t* seed) {
  Direction u = rotate_angle({uvw}, mu, phi, seed);
  uvw[0] = u.x;
  uvw[1] = u.y;
  uvw[2] = u.z;
}


Direction rotate_angle(Direction u, double mu, const double* phi, uint64_t* seed)
{
  // Sample azimuthal angle in [0,2pi) if none provided
  double phi_;
  if (phi != nullptr) {
    phi_ = (*phi);
  } else {
    phi_ = 2.0*PI*prn(seed);
  }

  // Precompute factors to save flops
  double sinphi = std::sin(phi_);
  double cosphi = std::cos(phi_);
  double a = std::sqrt(std::fmax(0., 1. - mu*mu));
  double b = std::sqrt(std::fmax(0., 1. - u.z*u.z));

  // Need to treat special case where sqrt(1 - w**2) is close to zero by
  // expanding about the v component rather than the w component
  if (b > 1e-10) {
    return {mu*u.x + a*(u.x*u.z*cosphi - u.y*sinphi) / b,
            mu*u.y + a*(u.y*u.z*cosphi + u.x*sinphi) / b,
            mu*u.z - a*b*cosphi};
  } else {
    b = std::sqrt(1. - u.y*u.y);
    return {mu*u.x + a*(u.x*u.y*cosphi + u.z*sinphi) / b,
            mu*u.y - a*b*cosphi,
            mu*u.z + a*(u.y*u.z*cosphi - u.x*sinphi) / b};
    // TODO: use the following code to make PolarAzimuthal distributions match
    // spherical coordinate conventions. Remove the related fixup code in
    // PolarAzimuthal::sample.
    //return {mu*u.x + a*(-u.x*u.y*sinphi + u.z*cosphi) / b,
    //        mu*u.y + a*b*sinphi,
    //        mu*u.z - a*(u.y*u.z*sinphi + u.x*cosphi) / b};
  }
}


double maxwell_spectrum(double T, uint64_t* seed) {
  // Set the random numbers
  double r1 = prn(seed);
  double r2 = prn(seed);
  double r3 = prn(seed);

  // determine cosine of pi/2*r
  double c = std::cos(PI / 2. * r3);

  // Determine outgoing energy
  double E_out = -T * (std::log(r1) + std::log(r2) * c * c);

  return E_out;
}


double normal_variate(double mean, double standard_deviation, uint64_t* seed) {
  // perhaps there should be a limit to the number of resamples
  while ( true ) {
    double v1 = 2 * prn(seed) - 1.;
    double v2 = 2 * prn(seed) - 1.;

    double r = std::pow(v1, 2) + std::pow(v2, 2);
    double r2 = std::pow(r, 2);
    if (r2 < 1) {
      double z = std::sqrt(-2.0 * std::log(r2)/r2);
      z *= (prn(seed) <= 0.5) ? v1 : v2;
      return mean + standard_deviation*z;
    }
  }
}

double muir_spectrum(double e0, double m_rat, double kt, uint64_t* seed) {
  // note sigma here is a factor of 2 shy of equation
  // 8 in https://permalink.lanl.gov/object/tr?what=info:lanl-repo/lareport/LA-05411-MS
  double sigma = std::sqrt(2.*e0*kt/m_rat);
  return normal_variate(e0, sigma, seed);
}


double watt_spectrum(double a, double b, uint64_t* seed) {
  double w = maxwell_spectrum(a, seed);
  double E_out = w + 0.25 * a * a * b + (2. * prn(seed) - 1.) * std::sqrt(a * a * b * w);

  return E_out;
}


void broaden_wmp_polynomials(double E, double dopp, int n, double factors[])
{
  // Factors is already pre-allocated
  double sqrtE = std::sqrt(E);
  double beta = sqrtE * dopp;
  double half_inv_dopp2 = 0.5 / (dopp * dopp);
  double quarter_inv_dopp4 = half_inv_dopp2 * half_inv_dopp2;

  double erf_beta;    // error function of beta
  double exp_m_beta2; // exp(-beta**2)
  if (beta > 6.0) {
    // Save time, ERF(6) is 1 to machine precision.
    // beta/sqrtpi*exp(-beta**2) is also approximately 1 machine epsilon.
    erf_beta = 1.;
    exp_m_beta2 = 0.;
  } else {
    erf_beta = std::erf(beta);
    exp_m_beta2 = std::exp(-beta * beta);
  }

  // Assume that, for sure, we'll use a second order (1/E, 1/V, const)
  // fit, and no less.

  factors[0] = erf_beta / E;
  factors[1] = 1. / sqrtE;
  factors[2] = factors[0] * (half_inv_dopp2 + E) + exp_m_beta2 /
       (beta * SQRT_PI);

  // Perform recursive broadening of high order components
  for (int i = 0; i < n - 3; i++) {
    double ip1_dbl = i + 1;
    if (i != 0) {
      factors[i + 3] = -factors[i - 1] * (ip1_dbl - 1.) * ip1_dbl *
           quarter_inv_dopp4 + factors[i + 1] *
           (E + (1. + 2. * ip1_dbl) * half_inv_dopp2);
    } else {
      // Although it's mathematically identical, factors[0] will contain
      // nothing, and we don't want to have to worry about memory.
      factors[i + 3] = factors[i + 1] *
           (E + (1. + 2. * ip1_dbl) * half_inv_dopp2);
    }
  }
}


void spline(int n, const double x[], const double y[], double z[])
{
  std::vector<double> c_new(n-1);

  // Set natural boundary conditions
  c_new[0] = 0.0;
  z[0] = 0.0;
  z[n-1] = 0.0;

  // Solve using tridiagonal matrix algorithm; first do forward sweep
  for (int i = 1; i < n - 1; i++) {
    double a = x[i] - x[i-1];
    double c = x[i+1] - x[i];
    double b = 2.0*(a + c);
    double d = 6.0*((y[i+1] - y[i])/c - (y[i] - y[i-1])/a);

    c_new[i] = c/(b - a*c_new[i-1]);
    z[i] = (d - a*z[i-1])/(b - a*c_new[i-1]);
  }

  // Back substitution
  for (int i = n - 2; i >= 0; i--) {
    z[i] = z[i] - c_new[i]*z[i+1];
  }
}


double spline_interpolate(int n, const double x[], const double y[],
  const double z[], double xint)
{
  // Find the lower bounding index in x of xint
  int i = n - 1;
  while (--i) {
    if (xint >= x[i]) break;
  }

  double h = x[i+1] - x[i];
  double r = xint - x[i];

  // Compute the coefficients
  double b = (y[i+1] - y[i])/h - (h/6.0)*(z[i+1] + 2.0*z[i]);
  double c = z[i]/2.0;
  double d = (z[i+1] - z[i])/(h*6.0);

  return y[i] + b*r + c*r*r + d*r*r*r;
}


double spline_integrate(int n, const double x[], const double y[],
  const double z[], double xa, double xb)
{
  // Find the lower bounding index in x of the lower limit of integration.
  int ia = n - 1;
  while (--ia) {
    if (xa >= x[ia]) break;
  }

  // Find the lower bounding index in x of the upper limit of integration.
  int ib = n - 1;
  while (--ib) {
    if (xb >= x[ib]) break;
  }

  // Evaluate the integral
  double s = 0.0;
  for (int i = ia; i <= ib; i++) {
    double h = x[i+1] - x[i];

    // Compute the coefficients
    double b = (y[i+1] - y[i])/h - (h/6.0)*(z[i+1] + 2.0*z[i]);
    double c = z[i]/2.0;
    double d = (z[i+1] - z[i])/(h*6.0);

    // Subtract the integral from x[ia] to xa
    if (i == ia) {
      double r = xa - x[ia];
      s = s - (y[i]*r + b/2.0*r*r + c/3.0*r*r*r + d/4.0*r*r*r*r);
    }

    // Integrate from x[ib] to xb in final interval
    if (i == ib) {
      h = xb - x[ib];
    }

    // Accumulate the integral
    s = s + y[i]*h + b/2.0*h*h + c/3.0*h*h*h + d/4.0*h*h*h*h;
  }

  return s;
}

std::complex<double> faddeeva(std::complex<double> z)
{
  // Technically, the value we want is given by the equation:
  // w(z) = I/pi * Integrate[Exp[-t^2]/(z-t), {t, -Infinity, Infinity}]
  // as shown in Equation 63 from Hwang, R. N. "A rigorous pole
  // representation of multilevel cross sections and its practical
  // applications." Nucl. Sci. Eng. 96.3 (1987): 192-209.
  //
  // The MIT Faddeeva function evaluates w(z) = exp(-z^2)erfc(-iz). These
  // two forms of the Faddeeva function are related by a transformation.
  //
  // If we call the integral form w_int, and the function form w_fun:
  // For imag(z) > 0, w_int(z) = w_fun(z)
  // For imag(z) < 0, w_int(z) = -conjg(w_fun(conjg(z)))

  // Note that Faddeeva::w will interpret zero as machine epsilon
  return z.imag() > 0.0 ? Faddeeva::w(z) :
    -std::conj(Faddeeva::w(std::conj(z)));
}

std::complex<double> w_derivative(std::complex<double> z, int order)
{
  using namespace std::complex_literals;
  switch (order) {
  case 0:
    return faddeeva(z);
  case 1:
    return -2.0*z*faddeeva(z) + 2.0i / SQRT_PI;
  default:
    return -2.0*z*w_derivative(z, order-1)
      - 2.0*(order-1)*w_derivative(z, order-2);
  }
}

} // namespace openmc
