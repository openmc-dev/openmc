
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

#endif // OPENMC_DAGMC_H

#endif

extern "C" bool dagmc_enabled();
