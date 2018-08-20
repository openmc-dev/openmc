#ifndef OPENMC_GEOMETRY_H
#define OPENMC_GEOMETRY_H

#include <cstdint>
#include <vector>

#include "particle.h"

namespace openmc {

extern "C" int openmc_root_universe;

extern std::vector<int64_t> overlap_check_count;

//==============================================================================
//! Check for overlapping cells at the particle's position.
//==============================================================================

extern "C" bool
check_cell_overlap(Particle* p);

//==============================================================================
//! Locate the particle in the geometry tree and set its geometry data fields.
//==============================================================================

extern "C" bool
find_cell(Particle* p, int n_search_cells, int* search_cells);

} // namespace openmc

#endif // OPENMC_GEOMETRY_H
