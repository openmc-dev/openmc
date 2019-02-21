#include "openmc/bank.h"

#include "openmc/capi.h"
#include "openmc/error.h"

#include <cstdint>

// Explicit template instantiation definition
template class std::vector<openmc::Bank>;

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace simulation {

int64_t n_bank;

std::vector<Bank> source_bank;
std::vector<Bank> fission_bank;
#ifdef _OPENMP
std::vector<Bank> master_fission_bank;
#endif

} // namespace simulation

//==============================================================================
// Non-member functions
//==============================================================================

void free_memory_bank()
{
  simulation::source_bank.clear();
#pragma omp parallel
  {
    simulation::fission_bank.clear();
  }
#ifdef _OPENMP
  simulation::master_fission_bank.clear();
#endif
}

//==============================================================================
// C API
//==============================================================================

extern "C" int openmc_source_bank(Bank** ptr, int64_t* n)
{
  if (simulation::source_bank.size() == 0) {
    set_errmsg("Source bank has not been allocated.");
    return OPENMC_E_ALLOCATE;
  } else {
    *ptr = simulation::source_bank.data();
    *n = simulation::source_bank.size();
    return 0;
  }
}

extern "C" int openmc_fission_bank(Bank** ptr, int64_t* n)
{
  if (simulation::fission_bank.size() == 0) {
    set_errmsg("Fission bank has not been allocated.");
    return OPENMC_E_ALLOCATE;
  } else {
    *ptr = simulation::fission_bank.data();
    *n = simulation::fission_bank.size();
    return 0;
  }
}

//==============================================================================
// Fortran compatibility
//==============================================================================

extern "C" int fission_bank_delayed_group(int64_t i) {
  return simulation::fission_bank[i-1].delayed_group;
}

extern "C" double fission_bank_E(int64_t i) {
  return simulation::fission_bank[i-1].E;
}

extern "C" double fission_bank_wgt(int64_t i) {
  return simulation::fission_bank[i-1].wgt;
}

extern "C" void source_bank_xyz(int64_t i, double* xyz)
{
  xyz[0] = simulation::source_bank[i-1].xyz[0];
  xyz[1] = simulation::source_bank[i-1].xyz[1];
  xyz[2] = simulation::source_bank[i-1].xyz[2];
}

extern "C" double source_bank_E(int64_t i)
{
  return simulation::source_bank[i-1].E;
}

extern "C" double source_bank_wgt(int64_t i)
{
  return simulation::source_bank[i-1].wgt;
}

extern "C" void source_bank_set_wgt(int64_t i, double wgt)
{
  simulation::source_bank[i-1].wgt = wgt;
}

} // namespace openmc
