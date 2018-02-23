#ifndef CAD_H
#define CAD_H

#include "DagMC.hpp"
#include "cell.h"
#include "surface.h"

extern "C" void load_cad_geometry_c();
extern "C" void dealloc_cad_c();

#endif // CAD_H
