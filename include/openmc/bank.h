#ifndef OPENMC_BANK_H
#define OPENMC_BANK_H

#include <cstdint>

#include "openmc/particle.h"
#include "openmc/position.h"
#include "openmc/shared_array.h"
#include "openmc/vector.h"

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace simulation {

extern vector<SourceSite> source_bank;

extern SharedArray<SourceSite> surf_source_bank;

extern SharedArray<SourceSite> fission_bank;

extern vector<int64_t> progeny_per_particle;

} // namespace simulation

//==============================================================================
// Non-member functions
//==============================================================================

void sort_fission_bank();

void free_memory_bank();

void init_fission_bank(int64_t max);

} // namespace openmc

#endif // OPENMC_BANK_H
