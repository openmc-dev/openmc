
#ifdef DAGMC

#ifndef DAGMC_H
#define DAGMC_H

#include "DagMC.hpp"
#include "openmc/cell.h"
#include "openmc/surface.h"

extern moab::DagMC* DAG;

extern "C" void load_dagmc_geometry_c();
extern "C" void free_memory_dagmc_c();

#endif // DAGMC_H

#endif
