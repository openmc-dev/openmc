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
extern std::vector<Particle::Bank> fission_bank;

extern Particle::Bank* shared_fission_bank;
extern int shared_fission_bank_length;
extern int shared_fission_bank_max;

} // namespace simulation

//==============================================================================
// Non-member functions
//==============================================================================

void free_memory_bank();

void init_shared_fission_bank(int max);

} // namespace openmc

#endif // OPENMC_BANK_H
