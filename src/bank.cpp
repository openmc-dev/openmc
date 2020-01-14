#include "openmc/bank.h"
#include "openmc/capi.h"
#include "openmc/error.h"

#include <cstdint>


namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace simulation {

std::vector<Particle::Bank> source_bank;

Particle::Bank* fission_bank {nullptr};
int64_t fission_bank_length {0};
int64_t fission_bank_max;

} // namespace simulation

//==============================================================================
// Non-member functions
//==============================================================================

void free_memory_bank()
{
  simulation::source_bank.clear();
  if( simulation::fission_bank != nullptr )
    delete[] simulation::fission_bank;
  simulation::fission_bank = nullptr;
  simulation::fission_bank_length = 0;
}

void init_fission_bank(int64_t max)
{
  simulation::fission_bank_max = max;
  simulation::fission_bank = new Particle::Bank[max];
  simulation::fission_bank_length = 0;
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
  if (simulation::fission_bank_length == 0) {
    set_errmsg("Fission bank has not been allocated.");
    return OPENMC_E_ALLOCATE;
  } else {
    *ptr = simulation::fission_bank;
    *n =   simulation::fission_bank_length;
    return 0;
  }
}

} // namespace openmc
