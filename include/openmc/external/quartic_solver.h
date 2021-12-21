#ifndef OPENMC_EXTERNAL_QUARTIC_SOLVER_H
#define OPENMC_EXTERNAL_QUARTIC_SOLVER_H

#include <complex>

namespace oqs {
void quartic_solver(double coeff[5], std::complex<double> roots[4]);
} // end namespace oqs

#endif // OPENMC_EXTERNAL_QUARTIC_SOLVER_H
