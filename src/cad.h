
#ifdef CAD

#ifndef CAD_H
#define CAD_H

#include "DagMC.hpp"
#include "cell.h"
#include "surface.h"

extern moab::DagMC* DAGMC;

extern "C" void load_cad_geometry_c();
extern "C" void free_memory_cad_c();

#endif // CAD_H

#endif
