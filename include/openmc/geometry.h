#ifndef OPENMC_GEOMETRY_H
#define OPENMC_GEOMETRY_H

#include <cstdint>
#include <vector>

#include "openmc/particle.h"


namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace model {

extern "C" int root_universe;

extern std::vector<int64_t> overlap_check_count;

} // namespace model

//==============================================================================
//! Check for overlapping cells at a particle's position.
//==============================================================================

extern "C" bool
check_cell_overlap(Particle* p);

//==============================================================================
//! Locate a particle in the geometry tree and set its geometry data fields.
//!
//! \param p A particle to be located.  This function will populate the
//!   geometry-dependent data fields of the particle.
//! \param search_surf A surface that the particle is expected to be on.  This
//!   value should be the signed, 1-based index of a surface.  If positive, the
//!   cells on the positive half-space of the surface will be searched.  If
//!   negative, the negative half-space will be searched.
//! \return True if the particle's location could be found and ascribed to a
//!   valid geometry coordinate stack.
//==============================================================================

extern "C" bool
find_cell(Particle* p, int search_surf);

//==============================================================================
//! Move a particle into a new lattice tile.
//==============================================================================

extern "C" void
cross_lattice(Particle* p, int lattice_translation[3]);

//==============================================================================
//! Find the next boundary a particle will intersect.
//==============================================================================

extern "C" void
distance_to_boundary(Particle* p, double* dist, int* surface_crossed,
                     int lattice_translation[3], int* next_level);

} // namespace openmc

#endif // OPENMC_GEOMETRY_H
