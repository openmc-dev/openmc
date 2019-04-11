#include "openmc/bank.h"

#include "openmc/capi.h"
#include "openmc/error.h"

#include <cstdint>

// Explicit template instantiation definition
template class std::vector<openmc::Particle::Bank>;

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace simulation {

std::vector<Particle::Bank> source_bank;
std::vector<Particle::Bank> fission_bank;
std::vector<Particle::Bank> secondary_bank;
#ifdef _OPENMP
std::vector<Particle::Bank> master_fission_bank;
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

extern "C" int openmc_source_bank(void** ptr, int64_t* n)
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

extern "C" int openmc_fission_bank(void** ptr, int64_t* n)
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

} // namespace openmc
