#include "openmc/bank.h"
#include "openmc/capi.h"
#include "openmc/error.h"
#include "openmc/simulation.h"

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

// Each entry in this vector corresponds to the number of progeny produced
// this generation for the particle located at that index. This vector is
// used to efficiently sort the fission bank after each iteration.
std::vector<int64_t> progeny_per_particle;

} // namespace simulation

//==============================================================================
// Non-member functions
//==============================================================================

void free_memory_bank()
{
  simulation::source_bank.clear();
  if (simulation::fission_bank)
    delete[] simulation::fission_bank;
  simulation::fission_bank = nullptr;
  simulation::fission_bank_length = 0;
}

void init_fission_bank(int64_t max)
{
  simulation::fission_bank_max = max;
  simulation::fission_bank = new Particle::Bank[max];
  simulation::fission_bank_length = 0;
  simulation::progeny_per_particle.resize(simulation::work_per_rank);
}

void sort_fission_bank()
{
  if( simulation::progeny_per_particle.size() == 0 )
    return;

  // Perform exclusive scan summation to determine starting indices in fission
  // bank for each parent particle id
  int64_t tmp = simulation::progeny_per_particle[0];
  simulation::progeny_per_particle[0] = 0;
  for (int64_t i = 1; i < simulation::progeny_per_particle.size(); i++) {
    int64_t value = simulation::progeny_per_particle[i-1] + tmp;
    tmp = simulation::progeny_per_particle[i];
    simulation::progeny_per_particle[i] = value;
  }

  // Create temporary scratch vector for permutation
  std::vector<Particle::Bank> sorted_bank(simulation::fission_bank_length);

  // Use parent and progeny indices to sort fission bank
  for (int64_t i = 0; i < simulation::fission_bank_length; i++) {
    Particle::Bank& site = simulation::fission_bank[i];
    int64_t idx = simulation::progeny_per_particle[site.parent_id-1] + site.progeny_id;
    sorted_bank[idx] = site;
  }

  // Copy sorted bank into the fission bank
  for (int64_t i = 0; i < simulation::fission_bank_length; i++) {
    simulation::fission_bank[i] = sorted_bank[i];
  }
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
