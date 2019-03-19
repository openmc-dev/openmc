#ifndef OPENMC_BANK_H
#define OPENMC_BANK_H

#include <cstdint>
#include <vector>

#include "openmc/particle.h"
#include "openmc/position.h"

// Without an explicit instantiation of vector<Bank>, the Intel compiler
// will complain about the threadprivate directive on filter_matches. Note that
// this has to happen *outside* of the openmc namespace
extern template class std::vector<openmc::Particle::Bank>;

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace simulation {

extern std::vector<Particle::Bank> source_bank;
extern std::vector<Particle::Bank> fission_bank;
extern std::vector<Particle::Bank> secondary_bank;
#ifdef _OPENMP
extern std::vector<Particle::Bank> master_fission_bank;
#endif

#pragma omp threadprivate(fission_bank, secondary_bank)

} // namespace simulation

//==============================================================================
// Non-member functions
//==============================================================================

void free_memory_bank();

} // namespace openmc

#endif // OPENMC_BANK_H
