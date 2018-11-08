
#ifndef OPENMC_DAGMC_H
#define OPENMC_DAGMC_H

#ifdef DAGMC

#include "DagMC.hpp"
#include "openmc/cell.h"
#include "openmc/surface.h"

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace model {

extern moab::DagMC* DAG;

} // namespace model

//==============================================================================
// Non-member functions
//==============================================================================

extern "C" void load_dagmc_geometry();
extern "C" void free_memory_dagmc();

} // namespace openmc

#endif // DAGMC

#endif // OPENMC_DAGMC_H

#ifdef DAGMC
extern "C" constexpr bool dagmc_enabled = true;
#else
extern "C" constexpr bool dagmc_enabled = false;
#endif
