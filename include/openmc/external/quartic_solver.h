#ifndef OPENMC_EXTERNAL_QUARTIC_SOLVER_H
#define OPENMC_EXTERNAL_QUARTIC_SOLVER_H

#include <complex>

void oqs_quartic_solver(double coeff[5], std::complex<double> roots[4]);

#endif // OPENMC_EXTERNAL_QUARTIC_SOLVER_H
