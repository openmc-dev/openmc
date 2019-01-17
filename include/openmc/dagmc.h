
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

extern "C" const bool dagmc_enabled;

namespace model {

extern moab::DagMC* DAG;

} // namespace model

//==============================================================================
// Non-member functions
//==============================================================================

extern "C" void load_dagmc_geometry();
extern "C" void free_memory_dagmc();
extern "C" pugi::xml_document* read_uwuw_materials();
bool get_uwuw_materials_xml(std::string& s);
} // namespace openmc

#endif // DAGMC

#endif // OPENMC_DAGMC_H
