#include <complex.h>
#include <float.h>
#include <math.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#define Sqr(x) ((x) * (x))
#ifndef CMPLX
#define CMPLX(x, y) (x) + (y)*I
#endif
const double cubic_rescal_fact =
  3.488062113727083E+102; //= pow(DBL_MAX,1.0/3.0)/1.618034;
const double quart_rescal_fact =
  7.156344627944542E+76; // = pow(DBL_MAX,1.0/4.0)/1.618034;
const double macheps = 2.2204460492503131E-16; // DBL_EPSILON
double oqs_max2(double a, double b)
{
  if (a >= b)
    return a;
  else
    return b;
}
double oqs_max3(double a, double b, double c)
{
  double t;
  t = oqs_max2(a, b);
  return oqs_max2(t, c);
}

void oqs_solve_cubic_analytic_depressed_handle_inf(
  double b, double c, double* sol)
{
  /* find analytically the dominant root of a depressed cubic x^3+b*x+c
   * where coefficients b and c are large (see sec. 2.2 in the manuscript) */
  double Q, R, theta, A, B, QR, QRSQ, KK, sqrtQ, RQ;
  ;
  const double PI2 = M_PI / 2.0, TWOPI = 2.0 * M_PI;
  Q = -b / 3.0;
  R = 0.5 * c;
  if (R == 0) {
    if (b <= 0) {
      *sol = sqrt(-b);
    } else {
      *sol = 0;
    }
    return;
  }

  if (fabs(Q) < fabs(R)) {
    QR = Q / R;
    QRSQ = QR * QR;
    KK = 1.0 - Q * QRSQ;
  } else {
    RQ = R / Q;
    KK = copysign(1.0, Q) * (RQ * RQ / Q - 1.0);
  }

  if (KK < 0.0) {
    sqrtQ = sqrt(Q);
    theta = acos((R / fabs(Q)) / sqrtQ);
    if (theta < PI2)
      *sol = -2.0 * sqrtQ * cos(theta / 3.0);
    else
      *sol = -2.0 * sqrtQ * cos((theta + TWOPI) / 3.0);
  } else {
    if (fabs(Q) < fabs(R))
      A = -copysign(1.0, R) * cbrt(fabs(R) * (1.0 + sqrt(KK)));
    else {
      A =
        -copysign(1.0, R) * cbrt(fabs(R) + sqrt(fabs(Q)) * fabs(Q) * sqrt(KK));
    }
    if (A == 0.0)
      B = 0.0;
    else
      B = Q / A;
    *sol = A + B;
  }
}
void oqs_solve_cubic_analytic_depressed(double b, double c, double* sol)
{
  /* find analytically the dominant root of a depressed cubic x^3+b*x+c
   * (see sec. 2.2 in the manuscript) */
  double Q, R, theta, Q3, R2, A, B, sqrtQ;
  Q = -b / 3.0;
  R = 0.5 * c;
  if (fabs(Q) > 1E102 || fabs(R) > 1E154) {
    oqs_solve_cubic_analytic_depressed_handle_inf(b, c, sol);
    return;
  }
  Q3 = Sqr(Q) * Q;
  R2 = Sqr(R);
  if (R2 < Q3) {
    theta = acos(R / sqrt(Q3));
    sqrtQ = -2.0 * sqrt(Q);
    if (theta < M_PI / 2)
      *sol = sqrtQ * cos(theta / 3.0);
    else
      *sol = sqrtQ * cos((theta + 2.0 * M_PI) / 3.0);
  } else {
    A = -copysign(1.0, R) * pow(fabs(R) + sqrt(R2 - Q3), 1.0 / 3.0);
    if (A == 0.0)
      B = 0.0;
    else
      B = Q / A;
    *sol = A + B; /* this is always largest root even if A=B */
  }
}
void oqs_calc_phi0(
  double a, double b, double c, double d, double* phi0, int scaled)
{
  /* find phi0 as the dominant root of the depressed and shifted cubic
   * in eq. (79) (see also the discussion in sec. 2.2 of the manuscript) */
  double rmax, g, h, gg, hh, aq, bq, cq, dq, s, diskr;
  double maxtt, xxx, gx, x, xold, f, fold, df, xsq;
  double ggss, hhss, dqss, aqs, bqs, cqs, rfact, rfactsq;
  int iter;
  diskr = 9 * a * a - 24 * b;
  /* eq. (87) */
  if (diskr > 0.0) {
    diskr = sqrt(diskr);
    if (a > 0.0)
      s = -2 * b / (3 * a + diskr);
    else
      s = -2 * b / (3 * a - diskr);
  } else {
    s = -a / 4;
  }
  /* eqs. (83) */
  aq = a + 4 * s;
  bq = b + 3 * s * (a + 2 * s);
  cq = c + s * (2 * b + s * (3 * a + 4 * s));
  dq = d + s * (c + s * (b + s * (a + s)));
  gg = bq * bq / 9;
  hh = aq * cq;

  g = hh - 4 * dq - 3 * gg;                                     /* eq. (85) */
  h = (8 * dq + hh - 2 * gg) * bq / 3 - cq * cq - dq * aq * aq; /* eq. (86) */
  oqs_solve_cubic_analytic_depressed(g, h, &rmax);
  if (isnan(rmax) || isinf(rmax)) {
    oqs_solve_cubic_analytic_depressed_handle_inf(g, h, &rmax);
    if ((isnan(rmax) || isinf(rmax)) && scaled) {
      // try harder: rescale also the depressed cubic if quartic has been
      // already rescaled
      rfact = cubic_rescal_fact;
      rfactsq = rfact * rfact;
      ggss = gg / rfactsq;
      hhss = hh / rfactsq;
      dqss = dq / rfactsq;
      aqs = aq / rfact;
      bqs = bq / rfact;
      cqs = cq / rfact;
      ggss = bqs * bqs / 9.0;
      hhss = aqs * cqs;
      g = hhss - 4.0 * dqss - 3.0 * ggss;
      h = (8.0 * dqss + hhss - 2.0 * ggss) * bqs / 3 - cqs * (cqs / rfact) -
          (dq / rfact) * aqs * aqs;
      oqs_solve_cubic_analytic_depressed(g, h, &rmax);
      if (isnan(rmax) || isinf(rmax)) {
        oqs_solve_cubic_analytic_depressed_handle_inf(g, h, &rmax);
      }
      rmax *= rfact;
    }
  }
  /* Newton-Raphson used to refine phi0 (see end of sec. 2.2 in the manuscript)
   */
  x = rmax;
  xsq = x * x;
  xxx = x * xsq;
  gx = g * x;
  f = x * (xsq + g) + h;
  if (fabs(xxx) > fabs(gx))
    maxtt = fabs(xxx);
  else
    maxtt = fabs(gx);
  if (fabs(h) > maxtt)
    maxtt = fabs(h);

  if (fabs(f) > macheps * maxtt) {
    for (iter = 0; iter < 8; iter++) {
      df = 3.0 * xsq + g;
      if (df == 0) {
        break;
      }
      xold = x;
      x += -f / df;
      fold = f;
      xsq = x * x;
      f = x * (xsq + g) + h;
      if (f == 0) {
        break;
      }

      if (fabs(f) >= fabs(fold)) {
        x = xold;
        break;
      }
    }
  }
  *phi0 = x;
}
double oqs_calc_err_ldlt(
  double b, double c, double d, double d2, double l1, double l2, double l3)
{
  /* Eqs. (29) and (30) in the manuscript */
  double sum;
  sum = (b == 0) ? fabs(d2 + l1 * l1 + 2.0 * l3)
                 : fabs(((d2 + l1 * l1 + 2.0 * l3) - b) / b);
  sum += (c == 0) ? fabs(2.0 * d2 * l2 + 2.0 * l1 * l3)
                  : fabs(((2.0 * d2 * l2 + 2.0 * l1 * l3) - c) / c);
  sum += (d == 0) ? fabs(d2 * l2 * l2 + l3 * l3)
                  : fabs(((d2 * l2 * l2 + l3 * l3) - d) / d);
  return sum;
}
double oqs_calc_err_abcd_cmplx(double a, double b, double c, double d,
  complex double aq, complex double bq, complex double cq, complex double dq)
{
  /* Eqs. (68) and (69) in the manuscript for complex alpha1 (aq), beta1 (bq),
   * alpha2 (cq) and beta2 (dq) */
  double sum;
  sum = (d == 0) ? cabs(bq * dq) : cabs((bq * dq - d) / d);
  sum +=
    (c == 0) ? cabs(bq * cq + aq * dq) : cabs(((bq * cq + aq * dq) - c) / c);
  sum +=
    (b == 0) ? cabs(bq + aq * cq + dq) : cabs(((bq + aq * cq + dq) - b) / b);
  sum += (a == 0) ? cabs(aq + cq) : cabs(((aq + cq) - a) / a);
  return sum;
}
double oqs_calc_err_abcd(double a, double b, double c, double d, double aq,
  double bq, double cq, double dq)
{
  /* Eqs. (68) and (69) in the manuscript for real alpha1 (aq), beta1 (bq),
   * alpha2 (cq) and beta2 (dq)*/
  double sum;
  sum = (d == 0) ? fabs(bq * dq) : fabs((bq * dq - d) / d);
  sum +=
    (c == 0) ? fabs(bq * cq + aq * dq) : fabs(((bq * cq + aq * dq) - c) / c);
  sum +=
    (b == 0) ? fabs(bq + aq * cq + dq) : fabs(((bq + aq * cq + dq) - b) / b);
  sum += (a == 0) ? fabs(aq + cq) : fabs(((aq + cq) - a) / a);
  return sum;
}
double oqs_calc_err_abc(
  double a, double b, double c, double aq, double bq, double cq, double dq)
{
  /* Eqs. (48)-(51) in the manuscript */
  double sum;
  sum =
    (c == 0) ? fabs(bq * cq + aq * dq) : fabs(((bq * cq + aq * dq) - c) / c);
  sum +=
    (b == 0) ? fabs(bq + aq * cq + dq) : fabs(((bq + aq * cq + dq) - b) / b);
  sum += (a == 0) ? fabs(aq + cq) : fabs(((aq + cq) - a) / a);
  return sum;
}
void oqs_NRabcd(double a, double b, double c, double d, double* AQ, double* BQ,
  double* CQ, double* DQ)
{
  /* Newton-Raphson described in sec. 2.3 of the manuscript for complex
   * coefficients a,b,c,d */
  int iter, k1, k2;
  double x02, errf, errfold, xold[4], x[4], dx[4], det, Jinv[4][4], fvec[4],
    vr[4];
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
  errf = 0;
  for (k1 = 0; k1 < 4; k1++) {
    errf += (vr[k1] == 0) ? fabs(fvec[k1]) : fabs(fvec[k1] / vr[k1]);
  }
  for (iter = 0; iter < 8; iter++) {
    x02 = x[0] - x[2];
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
    for (k1 = 0; k1 < 4; k1++) {
      dx[k1] = 0;
      for (k2 = 0; k2 < 4; k2++)
        dx[k1] += Jinv[k1][k2] * fvec[k2];
    }
    for (k1 = 0; k1 < 4; k1++)
      xold[k1] = x[k1];

    for (k1 = 0; k1 < 4; k1++) {
      x[k1] += -dx[k1] / det;
    }
    fvec[0] = x[1] * x[3] - d;
    fvec[1] = x[1] * x[2] + x[0] * x[3] - c;
    fvec[2] = x[1] + x[0] * x[2] + x[3] - b;
    fvec[3] = x[0] + x[2] - a;
    errfold = errf;
    errf = 0;
    for (k1 = 0; k1 < 4; k1++) {
      errf += (vr[k1] == 0) ? fabs(fvec[k1]) : fabs(fvec[k1] / vr[k1]);
    }
    if (errf == 0)
      break;
    if (errf >= errfold) {
      for (k1 = 0; k1 < 4; k1++)
        x[k1] = xold[k1];
      break;
    }
  }
  *AQ = x[0];
  *BQ = x[1];
  *CQ = x[2];
  *DQ = x[3];
}
void oqs_solve_quadratic(double a, double b, complex double roots[2])
{
  double div, sqrtd, diskr, zmax, zmin;
  diskr = a * a - 4 * b;
  if (diskr >= 0.0) {
    if (a >= 0.0)
      div = -a - sqrt(diskr);
    else
      div = -a + sqrt(diskr);

    zmax = div / 2;

    if (zmax == 0.0)
      zmin = 0.0;
    else
      zmin = b / zmax;

    roots[0] = CMPLX(zmax, 0.0);
    roots[1] = CMPLX(zmin, 0.0);
  } else {
    sqrtd = sqrt(-diskr);
    roots[0] = CMPLX(-a / 2, sqrtd / 2);
    roots[1] = CMPLX(-a / 2, -sqrtd / 2);
  }
}
void oqs_quartic_solver(double coeff[5], complex double roots[4])
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
  complex double acx1, bcx1, ccx1, dcx1, acx, bcx, ccx, dcx, cdiskr, zx1, zx2,
    zxmax, zxmin, qroots[2];
  double l2m[12], d2m[12], res[12], resmin, bl311, dml3l3, err0 = 0, err1 = 0,
                                                           aq1, bq1, cq1, dq1;
  double a, b, c, d, phi0, aq, bq, cq, dq, d2, d3, l1, l2, l3, errmin, errv[3],
    aqv[3], cqv[3], gamma, del2;
  int realcase[2], whichcase, k1, k, kmin, nsol;
  double rfactsq, rfact = 1.0;

  if (coeff[4] == 0.0) {
    printf("That's not a quartic!\n");
    return;
  }
  a = coeff[3] / coeff[4];
  b = coeff[2] / coeff[4];
  c = coeff[1] / coeff[4];
  d = coeff[0] / coeff[4];
  oqs_calc_phi0(a, b, c, d, &phi0, 0);

  // simple polynomial rescaling
  if (isnan(phi0) || isinf(phi0)) {
    rfact = quart_rescal_fact;
    a /= rfact;
    rfactsq = rfact * rfact;
    b /= rfactsq;
    c /= rfactsq * rfact;
    d /= rfactsq * rfactsq;
    oqs_calc_phi0(a, b, c, d, &phi0, 1);
  }
  l1 = a / 2;            /* eq. (16) */
  l3 = b / 6 + phi0 / 2; /* eq. (18) */
  del2 = c - a * l3;     /* defined just after eq. (27) */
  nsol = 0;
  bl311 = 2. * b / 3. - phi0 - l1 * l1; /* This is d2 as defined in eq. (20)*/
  dml3l3 = d - l3 * l3; /* dml3l3 is d3 as defined in eq. (15) with d2=0 */

  /* Three possible solutions for d2 and l2 (see eqs. (28) and discussion which
   * follows) */
  if (bl311 != 0.0) {
    d2m[nsol] = bl311;
    l2m[nsol] = del2 / (2.0 * d2m[nsol]);
    res[nsol] = oqs_calc_err_ldlt(b, c, d, d2m[nsol], l1, l2m[nsol], l3);
    nsol++;
  }
  if (del2 != 0) {
    l2m[nsol] = 2 * dml3l3 / del2;
    if (l2m[nsol] != 0) {
      d2m[nsol] = del2 / (2 * l2m[nsol]);
      res[nsol] = oqs_calc_err_ldlt(b, c, d, d2m[nsol], l1, l2m[nsol], l3);
      nsol++;
    }

    d2m[nsol] = bl311;
    l2m[nsol] = 2.0 * dml3l3 / del2;
    res[nsol] = oqs_calc_err_ldlt(b, c, d, d2m[nsol], l1, l2m[nsol], l3);
    nsol++;
  }

  if (nsol == 0) {
    l2 = d2 = 0.0;
  } else {
    /* we select the (d2,l2) pair which minimizes errors */
    for (k1 = 0; k1 < nsol; k1++) {
      if (k1 == 0 || res[k1] < resmin) {
        resmin = res[k1];
        kmin = k1;
      }
    }
    d2 = d2m[kmin];
    l2 = l2m[kmin];
  }
  whichcase = 0;
  if (d2 < 0.0) {
    /* Case I eqs. (37)-(40) */
    gamma = sqrt(-d2);
    aq = l1 + gamma;
    bq = l3 + gamma * l2;

    cq = l1 - gamma;
    dq = l3 - gamma * l2;
    if (fabs(dq) < fabs(bq))
      dq = d / bq;
    else if (fabs(dq) > fabs(bq))
      bq = d / dq;
    if (fabs(aq) < fabs(cq)) {
      nsol = 0;
      if (dq != 0) {
        aqv[nsol] = (c - bq * cq) / dq; /* see eqs. (47) */
        errv[nsol] = oqs_calc_err_abc(a, b, c, aqv[nsol], bq, cq, dq);
        nsol++;
      }
      if (cq != 0) {
        aqv[nsol] = (b - dq - bq) / cq; /* see eqs. (47) */
        errv[nsol] = oqs_calc_err_abc(a, b, c, aqv[nsol], bq, cq, dq);
        nsol++;
      }
      aqv[nsol] = a - cq; /* see eqs. (47) */
      errv[nsol] = oqs_calc_err_abc(a, b, c, aqv[nsol], bq, cq, dq);
      nsol++;
      /* we select the value of aq (i.e. alpha1 in the manuscript) which
       * minimizes errors */
      for (k = 0; k < nsol; k++) {
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
        errv[nsol] = oqs_calc_err_abc(a, b, c, aq, bq, cqv[nsol], dq);
        nsol++;
      }
      if (aq != 0) {
        cqv[nsol] = (b - bq - dq) / aq; /* see eqs. (53) */
        errv[nsol] = oqs_calc_err_abc(a, b, c, aq, bq, cqv[nsol], dq);
        nsol++;
      }
      cqv[nsol] = a - aq; /* see eqs. (53) */
      errv[nsol] = oqs_calc_err_abc(a, b, c, aq, bq, cqv[nsol], dq);
      nsol++;
      /* we select the value of cq (i.e. alpha2 in the manuscript) which
       * minimizes errors */
      for (k = 0; k < nsol; k++) {
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
    gamma = sqrt(d2);
    acx = CMPLX(l1, gamma);
    bcx = CMPLX(l3, gamma * l2);
    ccx = conj(acx);
    dcx = conj(bcx);
    realcase[0] = 0;
  } else
    realcase[0] = -1; // d2=0
  /* Case III: d2 is 0 or approximately 0 (in this case check which solution is
   * better) */
  if (realcase[0] == -1 || (fabs(d2) <= macheps * oqs_max3(fabs(2. * b / 3.),
                                                    fabs(phi0), l1 * l1))) {
    d3 = d - l3 * l3;
    if (realcase[0] == 1)
      err0 = oqs_calc_err_abcd(a, b, c, d, aq, bq, cq, dq);
    else if (realcase[0] == 0)
      err0 = oqs_calc_err_abcd_cmplx(a, b, c, d, acx, bcx, ccx, dcx);
    if (d3 <= 0) {
      realcase[1] = 1;
      aq1 = l1;
      bq1 = l3 + sqrt(-d3);
      cq1 = l1;
      dq1 = l3 - sqrt(-d3);
      if (fabs(dq1) < fabs(bq1))
        dq1 = d / bq1;
      else if (fabs(dq1) > fabs(bq1))
        bq1 = d / dq1;
      err1 = oqs_calc_err_abcd(a, b, c, d, aq1, bq1, cq1, dq1); /* eq. (68) */
    } else                                                      /* complex */
    {
      realcase[1] = 0;
      acx1 = l1;
      bcx1 = l3 + I * sqrt(d3);
      ccx1 = l1;
      dcx1 = conj(bcx1);
      err1 = oqs_calc_err_abcd_cmplx(a, b, c, d, acx1, bcx1, ccx1, dcx1);
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
    oqs_NRabcd(a, b, c, d, &aq, &bq, &cq, &dq);
    /* finally calculate the roots as roots of p1(x) and p2(x) (see end of
     * sec. 2.1) */
    oqs_solve_quadratic(aq, bq, qroots);
    roots[0] = qroots[0];
    roots[1] = qroots[1];
    oqs_solve_quadratic(cq, dq, qroots);
    roots[2] = qroots[0];
    roots[3] = qroots[1];
  } else {
    /* complex coefficients of p1 and p2 */
    if (whichcase == 0) // d2!=0
    {
      cdiskr = acx * acx / 4 - bcx;
      /* calculate the roots as roots of p1(x) and p2(x) (see end of sec. 2.1)
       */
      zx1 = -acx / 2 + csqrt(cdiskr);
      zx2 = -acx / 2 - csqrt(cdiskr);
      if (cabs(zx1) > cabs(zx2))
        zxmax = zx1;
      else
        zxmax = zx2;
      zxmin = bcx / zxmax;
      roots[0] = zxmin;
      roots[1] = conj(zxmin);
      roots[2] = zxmax;
      roots[3] = conj(zxmax);
    } else // d2 ~ 0
    {
      /* never gets here! */
      cdiskr = csqrt(acx * acx - 4.0 * bcx);
      zx1 = -0.5 * (acx + cdiskr);
      zx2 = -0.5 * (acx - cdiskr);
      if (cabs(zx1) > cabs(zx2))
        zxmax = zx1;
      else
        zxmax = zx2;
      zxmin = bcx / zxmax;
      roots[0] = zxmax;
      roots[1] = zxmin;
      cdiskr = csqrt(ccx * ccx - 4.0 * dcx);
      zx1 = -0.5 * (ccx + cdiskr);
      zx2 = -0.5 * (ccx - cdiskr);
      if (cabs(zx1) > cabs(zx2))
        zxmax = zx1;
      else
        zxmax = zx2;
      zxmin = dcx / zxmax;
      roots[2] = zxmax;
      roots[3] = zxmin;
    }
  }
  if (rfact != 1.0) {
    for (k = 0; k < 4; k++)
      roots[k] *= rfact;
  }
}
