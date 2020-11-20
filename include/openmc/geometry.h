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

extern int root_universe;  //!< Index of root universe
extern "C" int n_coord_levels; //!< Number of CSG coordinate levels

extern vector<int64_t> overlap_check_count;

} // namespace model

namespace gpu {
extern __constant__ int root_universe;
}

//==============================================================================
//! Check two distances by coincidence tolerance
//==============================================================================

HD inline bool coincident(double d1, double d2)
{
  return std::abs(d1 - d2) < FP_COINCIDENT;
}

//==============================================================================
//! Check for overlapping cells at a particle's position.
//==============================================================================

bool check_cell_overlap(Particle& p, bool error=true);

//==============================================================================
//! Locate a particle in the geometry tree and set its geometry data fields.
//!
//! \param p A particle to be located.  This function will populate the
//!   geometry-dependent data fields of the particle.
//! \return True if the particle's location could be found and ascribed to a
//!   valid geometry coordinate stack.
//==============================================================================
HD bool exhaustive_find_cell(Particle& p);
HD bool neighbor_list_find_cell(Particle& p); // Only usable on surface crossings

//==============================================================================
//! Move a particle into a new lattice tile.
//==============================================================================

HD void cross_lattice(Particle& p, const BoundaryInfo& boundary);

//==============================================================================
//! Find the next boundary a particle will intersect.
//==============================================================================

HD BoundaryInfo distance_to_boundary(Particle& p);

} // namespace openmc

#endif // OPENMC_GEOMETRY_H
