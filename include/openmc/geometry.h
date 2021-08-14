#ifndef OPENMC_GEOMETRY_H
#define OPENMC_GEOMETRY_H

#include <cmath>
#include <cstdint>

#include "openmc/array.h"
#include "openmc/constants.h"
#include "openmc/vector.h"

namespace openmc {

class BoundaryInfo;
class Particle;

//==============================================================================
// Global variables
//==============================================================================

namespace model {

extern int root_universe;      //!< Index of root universe
extern "C" int n_coord_levels; //!< Number of CSG coordinate levels

extern vector<int64_t> overlap_check_count;

} // namespace model

//==============================================================================
//! Check two distances by coincidence tolerance
//==============================================================================

inline bool coincident(double d1, double d2)
{
  return std::abs(d1 - d2) < FP_COINCIDENT;
}

//==============================================================================
//! Check for overlapping cells at a particle's position.
//==============================================================================

bool check_cell_overlap(Particle& p, bool error = true);

//==============================================================================
//! Get the cell instance for a particle at the specified universe level
//!
//! \param p A particle for which to compute the instance using
//!   its coordinates
//! \param level The level (zero indexed) of the geometry where the instance
//! should be computed. \return The instance of the cell at the specified level.
//==============================================================================

int cell_instance_at_level(const Particle& p, int level);

//==============================================================================
//! Locate a particle in the geometry tree and set its geometry data fields.
//!
//! \param p A particle to be located.  This function will populate the
//!   geometry-dependent data fields of the particle.
//! \return True if the particle's location could be found and ascribed to a
//!   valid geometry coordinate stack.
//==============================================================================
bool exhaustive_find_cell(Particle& p);
bool neighbor_list_find_cell(Particle& p); // Only usable on surface crossings

//==============================================================================
//! Move a particle into a new lattice tile.
//==============================================================================

void cross_lattice(Particle& p, const BoundaryInfo& boundary);

//==============================================================================
//! Find the next boundary a particle will intersect.
//==============================================================================

BoundaryInfo distance_to_boundary(Particle& p);

} // namespace openmc

#endif // OPENMC_GEOMETRY_H
