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

extern int root_universe;  //!< Index of root universe
extern int n_coord_levels; //!< Number of CSG coordinate levels

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
//! \param use_neighbor_lists If true, neighbor lists will be used to accelerate
//!   the geometry search, but this only works if the cell attribute of the
//!   particle's lowest coordinate level is valid and meaningful.
//! \return True if the particle's location could be found and ascribed to a
//!   valid geometry coordinate stack.
//==============================================================================

extern "C" bool
find_cell(Particle* p, bool use_neighbor_lists);

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
