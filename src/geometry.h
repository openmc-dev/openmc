#ifndef OPENMC_GEOMETRY_H
#define OPENMC_GEOMETRY_H

#include <cstdint>
#include <vector>

#include "particle.h"

namespace openmc {

extern "C" int openmc_root_universe;

extern std::vector<int64_t> overlap_check_count;

//==============================================================================
//! Check for overlapping cells at a particle's position.
//==============================================================================

extern "C" bool
check_cell_overlap(Particle* p);

//==============================================================================
//! Locate a particle in the geometry tree and set its geometry data fields.
//==============================================================================

extern "C" bool
find_cell(Particle* p, int search_surf);

//==============================================================================
//! Move a particle into a new lattice tile.
//==============================================================================

extern "C" void
cross_lattice(Particle* p, int lattice_translation[3]);

} // namespace openmc

#endif // OPENMC_GEOMETRY_H
