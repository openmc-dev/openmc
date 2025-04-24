#include <algorithm>
#define _USE_MATH_DEFINES // to make M_PI declared in Intel and MSVC compilers
#include <cmath>
#include <complex>
#include <cstdlib>
#include <limits>

namespace oqs {

// TODO: replace with <numbers> when we go for C++20
constexpr double PI {3.141592653589793238462643383279502884L};

// pow(DBL_MAX,1.0/3.0)/1.618034;
constexpr double CUBIC_RESCAL_FACT = 3.488062113727083e+102;
// pow(DBL_MAX,1.0/4.0)/1.618034;
constexpr double QUART_RESCAL_FACT = 7.156344627944542e+76;
constexpr double MACHEPS = std::numeric_limits<double>::epsilon();

double solve_cubic_analytic_depressed_handle_inf(double b, double c)
{
  /* find analytically the dominant root of a depressed cubic x^3+b*x+c
   * where coefficients b and c are large (see sec. 2.2 in the manuscript) */
  double Q = -b / 3.0;
  double R = 0.5 * c;
  if (R == 0) {
    return (b <= 0) ? std::sqrt(-b) : 0;
  }

  double KK;
  if (std::abs(Q) < std::abs(R)) {
    double QR = Q / R;
    double QRSQ = QR * QR;
    KK = 1.0 - Q * QRSQ;
  } else {
    double RQ = R / Q;
    KK = std::copysign(1.0, Q) * (RQ * RQ / Q - 1.0);
  }

  if (KK < 0.0) {
    double sqrtQ = std::sqrt(Q);
    double theta = std::acos((R / std::abs(Q)) / sqrtQ);
    if (2.0 * theta < PI)
      return -2.0 * sqrtQ * std::cos(theta / 3.0);
    else
      return -2.0 * sqrtQ * std::cos((theta + 2.0 * PI) / 3.0);
  } else {
    double A;
    if (std::abs(Q) < std::abs(R))
      A = -std::copysign(1.0, R) * cbrt(std::abs(R) * (1.0 + std::sqrt(KK)));
    else {
      A = -std::copysign(1.0, R) *
          cbrt(
            std::abs(R) + std::sqrt(std::abs(Q)) * std::abs(Q) * std::sqrt(KK));
    }
    double B = (A == 0.0) ? 0.0 : Q / A;
    return A + B;
  }
}

double solve_cubic_analytic_depressed(double b, double c)
{
  /* find analytically the dominant root of a depressed cubic x^3+b*x+c
   * (see sec. 2.2 in the manuscript) */
  double Q = -b / 3.0;
  double R = 0.5 * c;
  if (std::abs(Q) > 1e102 || std::abs(R) > 1e154) {
    return oqs::solve_cubic_analytic_depressed_handle_inf(b, c);
  }
  double Q3 = Q * Q * Q;
  double R2 = R * R;
  if (R2 < Q3) {
    double theta = std::acos(R / std::sqrt(Q3));
    double sqrtQ = -2.0 * std::sqrt(Q);
    if (2.0 * theta < PI)
      return sqrtQ * std::cos(theta / 3.0);
    else
      return sqrtQ * std::cos((theta + 2.0 * PI) / 3.0);
  } else {
    double A = -std::copysign(1.0, R) *
               std::pow(std::abs(R) + std::sqrt(R2 - Q3), 1.0 / 3.0);
    double B = (A == 0.0) ? 0.0 : Q / A;
    return A + B; /* this is always largest root even if A=B */
  }
}

double calc_phi0(double a, double b, double c, double d, int scaled)
{
  /* find phi0 as the dominant root of the depressed and shifted cubic
   * in eq. (79) (see also the discussion in sec. 2.2 of the manuscript) */
  double diskr = 9 * a * a - 24 * b;
  /* eq. (87) */
  double s;
  if (diskr > 0.0) {
    diskr = std::sqrt(diskr);
    s = -2 * b / (3 * a + std::copysign(diskr, a));
  } else {
    s = -a / 4;
  }
  /* eqs. (83) */
  double aq = a + 4 * s;
  double bq = b + 3 * s * (a + 2 * s);
  double cq = c + s * (2 * b + s * (3 * a + 4 * s));
  double dq = d + s * (c + s * (b + s * (a + s)));
  double gg = bq * bq / 9;
  double hh = aq * cq;

  double g = hh - 4 * dq - 3 * gg; /* eq. (85) */
  double h =
    (8 * dq + hh - 2 * gg) * bq / 3 - cq * cq - dq * aq * aq; /* eq. (86) */
  double rmax = oqs::solve_cubic_analytic_depressed(g, h);
  if (std::isnan(rmax) || std::isinf(rmax)) {
    rmax = oqs::solve_cubic_analytic_depressed_handle_inf(g, h);
    if ((std::isnan(rmax) || std::isinf(rmax)) && scaled) {
      // try harder: rescale also the depressed cubic if quartic has been
      // already rescaled
      double rfact = CUBIC_RESCAL_FACT;
      double rfactsq = rfact * rfact;
      double ggss = gg / rfactsq;
      double hhss = hh / rfactsq;
      double dqss = dq / rfactsq;
      double aqs = aq / rfact;
      double bqs = bq / rfact;
      double cqs = cq / rfact;
      ggss = bqs * bqs / 9.0;
      hhss = aqs * cqs;
      g = hhss - 4.0 * dqss - 3.0 * ggss;
      h = (8.0 * dqss + hhss - 2.0 * ggss) * bqs / 3 - cqs * (cqs / rfact) -
          (dq / rfact) * aqs * aqs;
      rmax = oqs::solve_cubic_analytic_depressed(g, h);
      if (std::isnan(rmax) || std::isinf(rmax)) {
        rmax = oqs::solve_cubic_analytic_depressed_handle_inf(g, h);
      }
      rmax *= rfact;
    }
  }
  /* Newton-Raphson used to refine phi0 (see end of sec. 2.2 in the manuscript)
   */
  double x = rmax;
  double xsq = x * x;
  double xxx = x * xsq;
  double gx = g * x;
  double f = x * (xsq + g) + h;
  double maxtt = std::max(std::abs(xxx), std::abs(gx));
  if (std::abs(h) > maxtt)
    maxtt = std::abs(h);

  if (std::abs(f) > MACHEPS * maxtt) {
    for (int iter = 0; iter < 8; iter++) {
      double df = 3.0 * xsq + g;
      if (df == 0) {
        break;
      }
      double xold = x;
      x += -f / df;
      double fold = f;
      xsq = x * x;
      f = x * (xsq + g) + h;
      if (f == 0) {
        break;
      }

      if (std::abs(f) >= std::abs(fold)) {
        x = xold;
        break;
      }
    }
  }
  return x;
}

double calc_err_ldlt(
  double b, double c, double d, double d2, double l1, double l2, double l3)
{
  /* Eqs. (29) and (30) in the manuscript */
  double sum = (b == 0) ? std::abs(d2 + l1 * l1 + 2.0 * l3)
                        : std::abs(((d2 + l1 * l1 + 2.0 * l3) - b) / b);
  sum += (c == 0) ? std::abs(2.0 * d2 * l2 + 2.0 * l1 * l3)
                  : std::abs(((2.0 * d2 * l2 + 2.0 * l1 * l3) - c) / c);
  sum += (d == 0) ? std::abs(d2 * l2 * l2 + l3 * l3)
                  : std::abs(((d2 * l2 * l2 + l3 * l3) - d) / d);
  return sum;
}

double calc_err_abcd_cmplx(double a, double b, double c, double d,
  std::complex<double> aq, std::complex<double> bq, std::complex<double> cq,
  std::complex<double> dq)
{
  /* Eqs. (68) and (69) in the manuscript for complex alpha1 (aq), beta1 (bq),
   * alpha2 (cq) and beta2 (dq) */
  double sum = (d == 0) ? std::abs(bq * dq) : std::abs((bq * dq - d) / d);
  sum += (c == 0) ? std::abs(bq * cq + aq * dq)
                  : std::abs(((bq * cq + aq * dq) - c) / c);
  sum += (b == 0) ? std::abs(bq + aq * cq + dq)
                  : std::abs(((bq + aq * cq + dq) - b) / b);
  sum += (a == 0) ? std::abs(aq + cq) : std::abs(((aq + cq) - a) / a);
  return sum;
}
double calc_err_abcd(double a, double b, double c, double d, double aq,
  double bq, double cq, double dq)
{
  /* Eqs. (68) and (69) in the manuscript for real alpha1 (aq), beta1 (bq),
   * alpha2 (cq) and beta2 (dq)*/
  double sum = (d == 0) ? std::abs(bq * dq) : std::abs((bq * dq - d) / d);
  sum += (c == 0) ? std::abs(bq * cq + aq * dq)
                  : std::abs(((bq * cq + aq * dq) - c) / c);
  sum += (b == 0) ? std::abs(bq + aq * cq + dq)
                  : std::abs(((bq + aq * cq + dq) - b) / b);
  sum += (a == 0) ? std::abs(aq + cq) : std::abs(((aq + cq) - a) / a);
  return sum;
}

double calc_err_abc(
  double a, double b, double c, double aq, double bq, double cq, double dq)
{
  /* Eqs. (48)-(51) in the manuscript */
  double sum = (c == 0) ? std::abs(bq * cq + aq * dq)
                        : std::abs(((bq * cq + aq * dq) - c) / c);
  sum += (b == 0) ? std::abs(bq + aq * cq + dq)
                  : std::abs(((bq + aq * cq + dq) - b) / b);
  sum += (a == 0) ? std::abs(aq + cq) : std::abs(((aq + cq) - a) / a);
  return sum;
}
void NRabcd(double a, double b, double c, double d, double* AQ, double* BQ,
  double* CQ, double* DQ)
{
  /* Newton-Raphson described in sec. 2.3 of the manuscript for complex
   * coefficients a,b,c,d */
  double xold[4], x[4], dx[4], det, Jinv[4][4], fvec[4], vr[4];
  x[0] = *AQ;
  x[1] = *BQ;
  x[2] = *CQ;
  x[3] = *DQ;
  vr[0] = d;
  vr[1] = c;
  vr[2] = b;
  vr[3] = a;
  fvec[0] = x[1] * x[3] - d;
  fvec[1] = x[1] * x[2] + x[0] * x[3] - c;
  fvec[2] = x[1] + x[0] * x[2] + x[3] - b;
  fvec[3] = x[0] + x[2] - a;
  double errf = 0;
  for (int k1 = 0; k1 < 4; k1++) {
    errf += (vr[k1] == 0) ? std::abs(fvec[k1]) : std::abs(fvec[k1] / vr[k1]);
  }
  for (int iter = 0; iter < 8; iter++) {
    double x02 = x[0] - x[2];
    det = x[1] * x[1] + x[1] * (-x[2] * x02 - 2.0 * x[3]) +
          x[3] * (x[0] * x02 + x[3]);
    if (det == 0.0)
      break;
    Jinv[0][0] = x02;
    Jinv[0][1] = x[3] - x[1];
    Jinv[0][2] = x[1] * x[2] - x[0] * x[3];
    Jinv[0][3] = -x[1] * Jinv[0][1] - x[0] * Jinv[0][2];
    Jinv[1][0] = x[0] * Jinv[0][0] + Jinv[0][1];
    Jinv[1][1] = -x[1] * Jinv[0][0];
    Jinv[1][2] = -x[1] * Jinv[0][1];
    Jinv[1][3] = -x[1] * Jinv[0][2];
    Jinv[2][0] = -Jinv[0][0];
    Jinv[2][1] = -Jinv[0][1];
    Jinv[2][2] = -Jinv[0][2];
    Jinv[2][3] = Jinv[0][2] * x[2] + Jinv[0][1] * x[3];
    Jinv[3][0] = -x[2] * Jinv[0][0] - Jinv[0][1];
    Jinv[3][1] = Jinv[0][0] * x[3];
    Jinv[3][2] = x[3] * Jinv[0][1];
    Jinv[3][3] = x[3] * Jinv[0][2];
    for (int k1 = 0; k1 < 4; k1++) {
      dx[k1] = 0;
      for (int k2 = 0; k2 < 4; k2++)
        dx[k1] += Jinv[k1][k2] * fvec[k2];
    }
    for (int k1 = 0; k1 < 4; k1++)
      xold[k1] = x[k1];

    for (int k1 = 0; k1 < 4; k1++) {
      x[k1] += -dx[k1] / det;
    }
    fvec[0] = x[1] * x[3] - d;
    fvec[1] = x[1] * x[2] + x[0] * x[3] - c;
    fvec[2] = x[1] + x[0] * x[2] + x[3] - b;
    fvec[3] = x[0] + x[2] - a;
    double errfold = errf;
    errf = 0;
    for (int k1 = 0; k1 < 4; k1++) {
      errf += (vr[k1] == 0) ? std::abs(fvec[k1]) : std::abs(fvec[k1] / vr[k1]);
    }
    if (errf == 0)
      break;
    if (errf >= errfold) {
      for (int k1 = 0; k1 < 4; k1++)
        x[k1] = xold[k1];
      break;
    }
  }
  *AQ = x[0];
  *BQ = x[1];
  *CQ = x[2];
  *DQ = x[3];
}

void solve_quadratic(double a, double b, std::complex<double> roots[2])
{
  double diskr = a * a - 4 * b;
  if (diskr >= 0.0) {
    double div = -a - std::copysign(std::sqrt(diskr), a);
    double zmax = div / 2;
    double zmin = (zmax == 0.0) ? 0.0 : b / zmax;

    roots[0] = std::complex<double>(zmax, 0.0);
    roots[1] = std::complex<double>(zmin, 0.0);
  } else {
    double sqrtd = std::sqrt(-diskr);
    roots[0] = std::complex<double>(-a / 2, sqrtd / 2);
    roots[1] = std::complex<double>(-a / 2, -sqrtd / 2);
  }
}

void quartic_solver(double coeff[5], std::complex<double> roots[4])
{
  /* USAGE:
   *
   * This routine calculates the roots of the quartic equation
   *
   * coeff[4]*x^4 + coeff[3]*x^3 + coeff[2]*x^2 + coeff[1]*x + coeff[0] = 0
   *
   * if coeff[4] != 0
   *
   * the four roots will be stored in the complex array roots[]
   *
   * */
  std::complex<double> acx, bcx, ccx, dcx;
  double l2m[12], d2m[12], res[12];
  double errv[3], aqv[3], cqv[3];
  int realcase[2];

  double a = coeff[3] / coeff[4];
  double b = coeff[2] / coeff[4];
  double c = coeff[1] / coeff[4];
  double d = coeff[0] / coeff[4];
  double phi0 = oqs::calc_phi0(a, b, c, d, 0);

  // simple polynomial rescaling
  double rfact = 1.0;
  if (std::isnan(phi0) || std::isinf(phi0)) {
    rfact = QUART_RESCAL_FACT;
    a /= rfact;
    double rfactsq = rfact * rfact;
    b /= rfactsq;
    c /= rfactsq * rfact;
    d /= rfactsq * rfactsq;
    phi0 = oqs::calc_phi0(a, b, c, d, 1);
  }
  double l1 = a / 2;            /* eq. (16) */
  double l3 = b / 6 + phi0 / 2; /* eq. (18) */
  double del2 = c - a * l3;     /* defined just after eq. (27) */
  int nsol = 0;
  double bl311 =
    2. * b / 3. - phi0 - l1 * l1; /* This is d2 as defined in eq. (20)*/
  double dml3l3 =
    d - l3 * l3; /* dml3l3 is d3 as defined in eq. (15) with d2=0 */

  /* Three possible solutions for d2 and l2 (see eqs. (28) and discussion which
   * follows) */
  if (bl311 != 0.0) {
    d2m[nsol] = bl311;
    l2m[nsol] = del2 / (2.0 * d2m[nsol]);
    res[nsol] = oqs::calc_err_ldlt(b, c, d, d2m[nsol], l1, l2m[nsol], l3);
    nsol++;
  }
  if (del2 != 0) {
    l2m[nsol] = 2 * dml3l3 / del2;
    if (l2m[nsol] != 0) {
      d2m[nsol] = del2 / (2 * l2m[nsol]);
      res[nsol] = oqs::calc_err_ldlt(b, c, d, d2m[nsol], l1, l2m[nsol], l3);
      nsol++;
    }

    d2m[nsol] = bl311;
    l2m[nsol] = 2.0 * dml3l3 / del2;
    res[nsol] = oqs::calc_err_ldlt(b, c, d, d2m[nsol], l1, l2m[nsol], l3);
    nsol++;
  }

  double l2, d2;
  if (nsol == 0) {
    l2 = d2 = 0.0;
  } else {
    /* we select the (d2,l2) pair which minimizes errors */
    double resmin;
    int kmin;
    for (int k1 = 0; k1 < nsol; k1++) {
      if (k1 == 0 || res[k1] < resmin) {
        resmin = res[k1];
        kmin = k1;
      }
    }
    d2 = d2m[kmin];
    l2 = l2m[kmin];
  }
  int whichcase = 0;
  double aq, bq, cq, dq;
  if (d2 < 0.0) {
    /* Case I eqs. (37)-(40) */
    double gamma = std::sqrt(-d2);
    aq = l1 + gamma;
    bq = l3 + gamma * l2;

    cq = l1 - gamma;
    dq = l3 - gamma * l2;
    if (std::abs(dq) < std::abs(bq))
      dq = d / bq;
    else if (std::abs(dq) > std::abs(bq))
      bq = d / dq;
    if (std::abs(aq) < std::abs(cq)) {
      nsol = 0;
      if (dq != 0) {
        aqv[nsol] = (c - bq * cq) / dq; /* see eqs. (47) */
        errv[nsol] = oqs::calc_err_abc(a, b, c, aqv[nsol], bq, cq, dq);
        nsol++;
      }
      if (cq != 0) {
        aqv[nsol] = (b - dq - bq) / cq; /* see eqs. (47) */
        errv[nsol] = oqs::calc_err_abc(a, b, c, aqv[nsol], bq, cq, dq);
        nsol++;
      }
      aqv[nsol] = a - cq; /* see eqs. (47) */
      errv[nsol] = oqs::calc_err_abc(a, b, c, aqv[nsol], bq, cq, dq);
      nsol++;
      /* we select the value of aq (i.e. alpha1 in the manuscript) which
       * minimizes errors */
      double errmin;
      int kmin;
      for (int k = 0; k < nsol; k++) {
        if (k == 0 || errv[k] < errmin) {
          kmin = k;
          errmin = errv[k];
        }
      }
      aq = aqv[kmin];
    } else {
      nsol = 0;
      if (bq != 0) {
        cqv[nsol] = (c - aq * dq) / bq; /* see eqs. (53) */
        errv[nsol] = oqs::calc_err_abc(a, b, c, aq, bq, cqv[nsol], dq);
        nsol++;
      }
      if (aq != 0) {
        cqv[nsol] = (b - bq - dq) / aq; /* see eqs. (53) */
        errv[nsol] = oqs::calc_err_abc(a, b, c, aq, bq, cqv[nsol], dq);
        nsol++;
      }
      cqv[nsol] = a - aq; /* see eqs. (53) */
      errv[nsol] = oqs::calc_err_abc(a, b, c, aq, bq, cqv[nsol], dq);
      nsol++;
      /* we select the value of cq (i.e. alpha2 in the manuscript) which
       * minimizes errors */
      double errmin;
      int kmin;
      for (int k = 0; k < nsol; k++) {
        if (k == 0 || errv[k] < errmin) {
          kmin = k;
          errmin = errv[k];
        }
      }
      cq = cqv[kmin];
    }
    realcase[0] = 1;
  } else if (d2 > 0) {
    /* Case II eqs. (53)-(56) */
    double gamma = std::sqrt(d2);
    acx = std::complex<double>(l1, gamma);
    bcx = std::complex<double>(l3, gamma * l2);
    ccx = std::conj(acx);
    dcx = std::conj(bcx);
    realcase[0] = 0;
  } else
    realcase[0] = -1; // d2=0
  /* Case III: d2 is 0 or approximately 0 (in this case check which solution is
   * better) */
  if (realcase[0] == -1 ||
      (std::abs(d2) <=
        MACHEPS * (std::abs(2. * b / 3.) + std::abs(phi0) + l1 * l1))) {
    double d3 = d - l3 * l3;
    double err0 = 0.0;
    if (realcase[0] == 1)
      err0 = oqs::calc_err_abcd(a, b, c, d, aq, bq, cq, dq);
    else if (realcase[0] == 0)
      err0 = oqs::calc_err_abcd_cmplx(a, b, c, d, acx, bcx, ccx, dcx);
    double aq1, bq1, cq1, dq1;
    std::complex<double> acx1, bcx1, ccx1, dcx1;
    double err1 = 0.0;
    if (d3 <= 0) {
      realcase[1] = 1;
      aq1 = l1;
      bq1 = l3 + std::sqrt(-d3);
      cq1 = l1;
      dq1 = l3 - std::sqrt(-d3);
      if (std::abs(dq1) < std::abs(bq1))
        dq1 = d / bq1;
      else if (std::abs(dq1) > std::abs(bq1))
        bq1 = d / dq1;
      err1 = oqs::calc_err_abcd(a, b, c, d, aq1, bq1, cq1, dq1); /* eq. (68) */
    } else {
      /* complex */
      realcase[1] = 0;
      acx1 = l1;
      bcx1 = l3 + std::complex<double>(0., std::sqrt(d3));
      ccx1 = l1;
      dcx1 = std::conj(bcx1);
      err1 = oqs::calc_err_abcd_cmplx(a, b, c, d, acx1, bcx1, ccx1, dcx1);
    }
    if (realcase[0] == -1 || err1 < err0) {
      whichcase = 1; // d2 = 0
      if (realcase[1] == 1) {
        aq = aq1;
        bq = bq1;
        cq = cq1;
        dq = dq1;
      } else {
        acx = acx1;
        bcx = bcx1;
        ccx = ccx1;
        dcx = dcx1;
      }
    }
  }
  if (realcase[whichcase] == 1) {
    /* if alpha1, beta1, alpha2 and beta2 are real first refine
     * the coefficient through a Newton-Raphson */
    oqs::NRabcd(a, b, c, d, &aq, &bq, &cq, &dq);
    /* finally calculate the roots as roots of p1(x) and p2(x) (see end of
     * sec. 2.1) */
    std::complex<double> qroots[2];
    oqs::solve_quadratic(aq, bq, qroots);
    roots[0] = qroots[0];
    roots[1] = qroots[1];
    oqs::solve_quadratic(cq, dq, qroots);
    roots[2] = qroots[0];
    roots[3] = qroots[1];
  } else {
    /* complex coefficients of p1 and p2 */
    if (whichcase == 0) { // d2!=0
      auto cdiskr = 0.25 * acx * acx - bcx;
      /* calculate the roots as roots of p1(x) and p2(x) (see end of sec. 2.1)
       */
      auto zx1 = -0.5 * acx + std::sqrt(cdiskr);
      auto zx2 = -0.5 * acx - std::sqrt(cdiskr);
      auto zxmax = (std::abs(zx1) > std::abs(zx2)) ? zx1 : zx2;
      auto zxmin = bcx / zxmax;
      roots[0] = zxmin;
      roots[1] = std::conj(zxmin);
      roots[2] = zxmax;
      roots[3] = std::conj(zxmax);
    } else { // d2 ~ 0
      /* never gets here! */
      auto cdiskr = std::sqrt(acx * acx - 4.0 * bcx);
      auto zx1 = -0.5 * (acx + cdiskr);
      auto zx2 = -0.5 * (acx - cdiskr);
      auto zxmax = (std::abs(zx1) > std::abs(zx2)) ? zx1 : zx2;
      auto zxmin = bcx / zxmax;
      roots[0] = zxmax;
      roots[1] = zxmin;
      cdiskr = std::sqrt(ccx * ccx - 4.0 * dcx);
      zx1 = -0.5 * (ccx + cdiskr);
      zx2 = -0.5 * (ccx - cdiskr);
      zxmax = (std::abs(zx1) > std::abs(zx2)) ? zx1 : zx2;
      zxmin = dcx / zxmax;
      roots[2] = zxmax;
      roots[3] = zxmin;
    }
  }
  if (rfact != 1.0) {
    for (int k = 0; k < 4; k++)
      roots[k] *= rfact;
  }
}

} // namespace oqs
