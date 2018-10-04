
#ifndef OPENMC_DAGMC_H
#define OPENMC_DAGMC_H

#ifdef DAGMC

#include "DagMC.hpp"
#include "openmc/cell.h"
#include "openmc/surface.h"

namespace openmc {

extern moab::DagMC* DAG;

extern "C" void load_dagmc_geometry();
extern "C" void free_memory_dagmc();

}

#endif // DAGMC

#endif // OPENMC_DAGMC_H

#ifdef DAGMC
extern "C" constexpr bool dagmc_enabled = true;
#else
extern "C" constexpr bool dagmc_enabled = false;
#endif
