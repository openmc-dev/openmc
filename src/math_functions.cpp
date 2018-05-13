#include "math_functions.h"

namespace openmc {


//==============================================================================
// NORMAL_PERCENTILE calculates the percentile of the standard normal
// distribution with a specified probability level
//==============================================================================

double __attribute__ ((const)) normal_percentile_c(const double p) {
  double z;
  double q;
  double r;
  const double p_low = 0.02425;
  const double a[6] = {-3.969683028665376e1, 2.209460984245205e2,
                       -2.759285104469687e2, 1.383577518672690e2,
                       -3.066479806614716e1, 2.506628277459239e0};
  const double b[5] = {-5.447609879822406e1, 1.615858368580409e2,
                       -1.556989798598866e2, 6.680131188771972e1,
                       -1.328068155288572e1};
  const double c[6] = {-7.784894002430293e-3, -3.223964580411365e-1,
                       -2.400758277161838, -2.549732539343734,
                       4.374664141464968, 2.938163982698783};
  const double d[4] = {7.784695709041462e-3, 3.224671290700398e-1,
                       2.445134137142996, 3.754408661907416};

  // The rational approximation used here is from an unpublished work at
  // http://home.online.no/~pjacklam/notes/invnorm/

  if (p < p_low) {
    // Rational approximation for lower region.

    q = std::sqrt(-2.0 * std::log(p));
    z = (((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) /
          ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1.0);

  } else if (p <= 1.0 - p_low) {
    // Rational approximation for central region

    q = p - 0.5;
    r = q * q;
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

//==============================================================================
// T_PERCENTILE calculates the percentile of the Student's t distribution with
// a specified probability level and number of degrees of freedom
//==============================================================================

double __attribute__ ((const)) t_percentile_c(const double p, const int df){
  double t;
  double n;
  double k;
  double z;
  double z2;

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

    n = static_cast<double>(df);
    k = 1. / (n - 2.);
    z = normal_percentile_c(p);
    z2 = z * z;
    t = std::sqrt(n * k) * (z + (z2 - 3.) * z * k / 4. + ((5. * z2 - 56.) * z2 +
         75.) * z * k * k / 96. + (((z2 - 27.) * 3. * z2 + 417.) * z2 - 315.) *
         z * k * k * k / 384.);
  }

  return t;
}

//==============================================================================
// CALC_PN calculates the n-th order Legendre polynomials at the value of x.
//==============================================================================

void calc_pn_c(const int n, const double x, double pnx[]) {
  int l;

  pnx[0] = 1.;
  if (n >= 1) {
    pnx[1] = x;
  }

  // Use recursion relation to build the higher orders
  for (l = 1; l < n; l ++) {
    pnx[l + 1] = (static_cast<double>(2 * l + 1) * x * pnx[l] -
                  static_cast<double>(l) * pnx[l - 1]) /
         (static_cast<double>(l + 1));
  }
}

//==============================================================================
// EVALUATE_LEGENDRE Find the value of f(x) given a set of Legendre coefficients
// and the value of x
//==============================================================================

double __attribute__ ((const)) evaluate_legendre_c(const int n,
                                                   const double data[],
                                                   const double x) {
  double val;
  int l;
  double* pnx = new double[n + 1];

  val = 0.0;
  calc_pn_c(n, x, pnx);
  for (l = 0; l <= n; l++) {
    val += (static_cast<double>(l) + 0.5) * data[l] * pnx[l];
  }
  delete[] pnx;
  return val;
}

//==============================================================================
// CALC_RN calculates the n-th order spherical harmonics for a given angle
// (in terms of (u,v,w)).  All Rn,m values are provided (where -n<=m<=n) for n
// all 0 <= n
//==============================================================================

void calc_rn_c(const int n, const double uvw[3], double rn[]){
  double phi;
  double w;
  double w2m1;
  int i;
  int l;

  // rn[] is assumed to have already been allocated to the correct size

  // Store the cosine of the polar angle and the azimuthal angle
  w = uvw[2];
  if (uvw[0] == 0.) {
    phi = 0.;
  } else {
    phi = std::atan2(uvw[1], uvw[0]);
  }

  // Store the shorthand of 1-w * w
  w2m1 = 1. - w * w;

  // Now evaluate the spherical harmonics function
  rn[0] = 1.;
  i = 0;
  for (l = 1; l <= n; l++) {
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

//==============================================================================
// CALC_ZN calculates the n-th order modified Zernike polynomial moment for a
// given angle (rho, theta) location in the unit disk. The normalization of the
// polynomials is such tha the integral of Z_pq*Z_pq over the unit disk is
// exactly pi
//==============================================================================

void calc_zn_c(const int n, const double rho, const double phi, double zn[]) {
  // This procedure uses the modified Kintner's method for calculating Zernike
  // polynomials as outlined in Chong, C. W., Raveendran, P., & Mukundan,
  // R. (2003). A comparative analysis of algorithms for fast computation of
  // Zernike moments. Pattern Recognition, 36(3), 731-742.
  double sin_phi;              // Cosine of phi
  double cos_phi;              // Sine of phi
  double sin_phi_vec[n + 1];   // Sin[n * phi]
  double cos_phi_vec[n + 1];   // Cos[n * phi]
  double zn_mat[n + 1][n + 1]; // Matrix forms of the coefficients which are
                               // easier to work with
  // Variables for R_m_n calculation
  double k1;
  double k2;
  double k3;
  double k4;
  // Loop counters
  int i;
  int p;
  int q;

  // n == radial degree
  // m == azimuthal frequency

  // ===========================================================================
  // Determine vector of sin(n*phi) and cos(n*phi). This takes advantage of the
  // following recurrence relations so that only a single sin/cos have to be
  // evaluated (http://mathworld.wolfram.com/Multiple-AngleFormulas.html)
  //
  // sin(nx) = 2 cos(x) sin((n-1)x) - sin((n-2)x)
  // cos(nx) = 2 cos(x) cos((n-1)x) - cos((n-2)x)

  sin_phi = std::sin(phi);
  cos_phi = std::cos(phi);

  sin_phi_vec[0] = 1.0;
  cos_phi_vec[0] = 1.0;

  sin_phi_vec[1] = 2.0 * cos_phi;
  cos_phi_vec[1] = cos_phi;

  for (i = 2; i <= n; i++) {
    sin_phi_vec[i] = 2. * cos_phi * sin_phi_vec[i - 1] - sin_phi_vec[i - 2];
    cos_phi_vec[i] = 2. * cos_phi * cos_phi_vec[i - 1] - cos_phi_vec[i - 2];
  }

  for (i = 0; i <= n; i++) {
    sin_phi_vec[i] *=  sin_phi;
  }

  // ===========================================================================
  // Calculate R_pq(rho)

  // Fill the main diagonal first (Eq 3.9 in Chong)
  for (p = 0; p <= n; p++) {
    zn_mat[p][p] = std::pow(rho, p);
  }

  // Fill the 2nd diagonal (Eq 3.10 in Chong)
  for (q = 0; q <= n - 2; q++) {
    zn_mat[q][q+2] = (q + 2) * zn_mat[q+2][q+2] - (q + 1) * zn_mat[q][q];
  }

  // Fill in the rest of the values using the original results (Eq. 3.8 in Chong)
  for (p = 4; p <= n; p++) {
    k2 = static_cast<double> (2 * p * (p - 1) * (p - 2));
    for (q = p - 4; q >= 0; q -= 2) {
      k1 = static_cast<double>((p + q) * (p - q) * (p - 2)) / 2.;
      k3 = static_cast<double>(-q * q * (p - 1) - p * (p - 1) * (p - 2));
      k4 = static_cast<double>(-p * (p + q - 2) * (p - q - 2)) / 2.;
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
  i = 0;
  for (p = 0; p <= n; p++) {
    for (q = -p; q <= p; q += 2) {
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

//==============================================================================
// ROTATE_ANGLE rotates direction std::cosines through a polar angle whose
// cosine is mu and through an azimuthal angle sampled uniformly. Note that
// this is done with direct sampling rather than rejection as is done in MCNP
// and SERPENT.
//==============================================================================

void rotate_angle_c(double uvw[3], const double mu, double* phi) {
  double phi_;   // azimuthal angle
  double sinphi; // std::sine of azimuthal angle
  double cosphi; // cosine of azimuthal angle
  double a;      // sqrt(1 - mu^2)
  double b;      // sqrt(1 - w^2)
  double u0;     // original std::cosine in x direction
  double v0;     // original std::cosine in y direction
  double w0;     // original std::cosine in z direction

  // Copy original directional std::cosines
  u0 = uvw[0];
  v0 = uvw[1];
  w0 = uvw[2];

  // Sample azimuthal angle in [0,2pi) if none provided
  if (phi != NULL) {
    phi_ = (*phi);
  } else {
    phi_ = 2. * PI * prn();
  }

  // Precompute factors to save flops
  sinphi = std::sin(phi_);
  cosphi = std::cos(phi_);
  a = std::sqrt(std::max(0., 1. - mu * mu));
  b = std::sqrt(std::max(0., 1. - w0 * w0));

  // Need to treat special case where sqrt(1 - w**2) is close to zero by
  // expanding about the v component rather than the w component
  if (b > 1e-10) {
    uvw[0] = mu * u0 + a * (u0 * w0 * cosphi - v0 * sinphi) / b;
    uvw[1] = mu * v0 + a * (v0 * w0 * cosphi + u0 * sinphi) / b;
    uvw[2] = mu * w0 - a * b * cosphi;
  } else {
    b = std::sqrt(1. - v0 * v0);
    uvw[0] = mu * u0 + a * (u0 * v0 * cosphi + w0 * sinphi) / b;
    uvw[1] = mu * v0 - a * b * cosphi;
    uvw[2] = mu * w0 + a * (v0 * w0 * cosphi - u0 * sinphi) / b;
  }
}

//==============================================================================
// MAXWELL_SPECTRUM samples an energy from the Maxwell fission distribution
// based on a direct sampling scheme. The probability distribution function for
// a Maxwellian is given as p(x) = 2/(T*sqrt(pi))*sqrt(x/T)*exp(-x/T).
// This PDF can be sampled using rule C64 in the Monte Carlo Sampler LA-9721-MS.
//==============================================================================

double maxwell_spectrum_c(const double T) {
  double E_out; // Sampled Energy

  double r1;
  double r2;
  double r3;  // random numbers
  double c;   // cosine of pi/2*r3

  r1 = prn();
  r2 = prn();
  r3 = prn();

  // determine cosine of pi/2*r
  c = std::cos(PI / 2. * r3);

  // determine outgoing energy
  E_out = -T * (std::log(r1) + std::log(r2) * c * c);

  return E_out;
}

//==============================================================================
// WATT_SPECTRUM samples the outgoing energy from a Watt energy-dependent
// fission spectrum. Although fitted parameters exist for many nuclides,
// generally the continuous tabular distributions (LAW 4) should be used in
// lieu of the Watt spectrum. This direct sampling scheme is an unpublished
// scheme based on the original Watt spectrum derivation (See F. Brown's
// MC lectures).
//==============================================================================

double watt_spectrum_c(const double a, const double b) {
  double E_out; // Sampled Energy
  double w;     // sampled from Maxwellian

  w = maxwell_spectrum_c(a);
  E_out = w + 0.25 * a * a * b + (2. * prn() - 1.) * std::sqrt(a * a * b * w);

  return E_out;
}

//==============================================================================
// FADDEEVA the Faddeeva function, using Stephen Johnson's implementation
//==============================================================================

// double complex __attribute__ ((const)) faddeeva_c(double complex z) {
//   double complex wv; // The resultant w(z) value
//   double relerr;      // Target relative error in the inner loop of MIT Faddeeva

//   // Technically, the value we want is given by the equation:
//   // w(z) = I/Pi * Integrate[Exp[-t^2]/(z-t), {t, -Infinity, Infinity}]
//   // as shown in Equation 63 from Hwang, R. N. "A rigorous pole
//   // representation of multilevel cross sections and its practical
//   // applications." Nuclear Science and Engineering 96.3 (1987): 192-209.
//   //
//   // The MIT Faddeeva function evaluates w(z) = exp(-z^2)erfc(-iz). These
//   // two forms of the Faddeeva function are related by a transformation.
//   //
//   // If we call the integral form w_int, and the function form w_fun:
//   // For imag(z) > 0, w_int(z) = w_fun(z)
//   // For imag(z) < 0, w_int(z) = -conjg(w_fun(conjg(z)))

//   // Note that faddeeva_w will interpret zero as machine epsilon

//   relerr = 0.;
//   if (cimag(z) > 0.) {
//       wv = Faddeeva_w(z, relerr);
//   } else {
//     wv = -conj(Faddeeva_w(conj(z), relerr));
//   }

//   return wv;
// }

// double complex w_derivative_c(double complex z, int order){
//   double complex wv; // The resultant w(z) value

//   const double complex twoi_sqrtpi = 2.0 / SQRT_PI * I;

//   switch(order) {
//     case 0:
//       wv = faddeeva_c(z);
//       break;
//     case 1:
//       wv = -2. * z * faddeeva_c(z) + twoi_sqrtpi;
//       break;
//     default:
//       wv = -2. * z * w_derivative_c(z, order - 1) - 2. * (order - 1) *
//            w_derivative_c(z, order - 2);
//   }

//   return wv;
// }

//==============================================================================
// BROADEN_WMP_POLYNOMIALS Doppler broadens the windowed multipole curvefit.
// The curvefit is a polynomial of the form a/E + b/sqrt(E) + c + d sqrt(E) ...
//==============================================================================

void broaden_wmp_polynomials_c(const double E, const double dopp, const int n,
                               double factors[]) {
  // Factors is already pre-allocated

  double sqrtE;               // sqrt(energy)
  double beta;                // sqrt(atomic weight ratio * E / kT)
  double half_inv_dopp2;      // 0.5 / dopp**2
  double quarter_inv_dopp4;   // 0.25 / dopp**4
  double erf_beta;            // error function of beta
  double exp_m_beta2;         // exp(-beta**2)
  int i;
  double ip1_dbl;

  sqrtE = std::sqrt(E);
  beta = sqrtE * dopp;
  half_inv_dopp2 = 0.5 / (dopp * dopp);
  quarter_inv_dopp4 = half_inv_dopp2 * half_inv_dopp2;

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
  for (i = 0; i < n - 3; i++) {
    ip1_dbl = static_cast<double>(i + 1);
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

} // namespace openmc
