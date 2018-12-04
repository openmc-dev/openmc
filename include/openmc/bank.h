#ifndef OPENMC_BANK_H
#define OPENMC_BANK_H

#include <cstdint>
#include <vector>

namespace openmc {

struct Bank {
  double wgt;
  double xyz[3];
  double uvw[3];
  double E;
  int delayed_group;
  int particle;
};

} // namespace openmc

// Without an explicit instantiation of vector<Bank>, the Intel compiler
// will complain about the threadprivate directive on filter_matches. Note that
// this has to happen *outside* of the openmc namespace
extern template class std::vector<openmc::Bank>;

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace simulation {

extern "C" int64_t n_bank;

extern std::vector<Bank> source_bank;
extern std::vector<Bank> fission_bank;
#ifdef _OPENMP
extern std::vector<Bank> master_fission_bank;
#endif

#pragma omp threadprivate(fission_bank, n_bank)

} // namespace simulation

} // namespace openmc

#endif // OPENMC_BANK_H
