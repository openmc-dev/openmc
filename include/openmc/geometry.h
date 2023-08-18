#ifndef OPENMC_GEOMETRY_H
#define OPENMC_GEOMETRY_H

#include <cmath>
#include <cstdint>

#include "openmc/array.h"
#include "openmc/constants.h"
#include "openmc/vector.h"

namespace openmc {

class BoundaryInfo;
class Geometron;

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

bool check_cell_overlap(Geometron& p, bool error = true);

//==============================================================================
//! Get the cell instance for a particle at the specified universe level
//!
//! \param p A particle for which to compute the instance using
//!   its coordinates
//! \param level The level (zero indexed) of the geometry where the instance
//! should be computed. \return The instance of the cell at the specified level.
//==============================================================================

int cell_instance_at_level(const Geometron& p, int level);

//==============================================================================
//! Locate a particle in the geometry tree and set its geometry data fields.
//!
//! \param p A particle to be located.  This function will populate the
//!   geometry-dependent data fields of the particle.
//! \return True if the particle's location could be found and ascribed to a
//!   valid geometry coordinate stack.
//==============================================================================
bool exhaustive_find_cell(Geometron& p, bool verbose = false);
bool neighbor_list_find_cell(
  Geometron& p, bool verbose = false); // Only usable on surface crossings

//==============================================================================
//! Move a particle into a new lattice tile.
//==============================================================================

void cross_lattice(
  Geometron& p, const BoundaryInfo& boundary, bool verbose = false);

//==============================================================================
//! Find the next boundary a particle will intersect.
//==============================================================================

BoundaryInfo distance_to_boundary(Geometron& p);

/* Geometry routines can potentially throw this type of exception
 */
class ParticleLost : public std::exception {
public:
  enum class Reason {
    negative_lattice_distance,   // should not have negative distance to lattice
    bad_boundary_crossing,       // could not locate after crossing a boundary
    no_universe_outside_lattice, // undefined behavior. All space should be
                                 // defined.
    no_dagmc_intersection
  } reason;

  // Extra data to be passed up for exception handling or for detailed
  // error reporting. For instance no_dagmc_intersection will give the
  // DAGMC cell ID.
  int id;

  ParticleLost(Reason reason_a, int id_a = 0) : reason(reason_a), id(id_a) {}

  // Handles uncaught exception messages, e.g. Geometron with no ID.
  // These will not show up in the usual particle tracking, as more
  // detailed error reporting is provided in the Particle class.
  const char* what()
  {
    switch (reason) {
    case Reason::negative_lattice_distance:
      return "Negative distance to a lattice";
    case Reason::bad_boundary_crossing:
      return "Crossed a boundary and lost particle";
    case Reason::no_universe_outside_lattice:
      return "Outside lattice but no outer region defined";
    case Reason::no_dagmc_intersection:
      return "No intersection found with DAGMC cell";
    }
  }
};

} // namespace openmc

#endif // OPENMC_GEOMETRY_H
