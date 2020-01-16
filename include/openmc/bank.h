#ifndef OPENMC_BANK_H
#define OPENMC_BANK_H

#include <cstdint>
#include <vector>

#include "openmc/particle.h"
#include "openmc/position.h"

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace simulation {

extern std::vector<Particle::Bank> source_bank;

extern Particle::Bank* fission_bank;
extern int64_t fission_bank_length;
extern int64_t fission_bank_max;
extern std::vector<int64_t> progeny_per_particle;

} // namespace simulation

//==============================================================================
// Non-member functions
//==============================================================================

void sort_fission_bank();

void free_memory_bank();

void init_fission_bank(int64_t max);

} // namespace openmc

#endif // OPENMC_BANK_H
