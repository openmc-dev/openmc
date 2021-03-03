
#ifndef OPENMC_DAGMC_H
#define OPENMC_DAGMC_H

namespace openmc {
extern "C" const bool dagmc_enabled;
}

#ifdef DAGMC

#include "DagMC.hpp"
#include "openmc/xml_interface.h"
#include "openmc/position.h"

namespace openmc {

namespace model {
  extern moab::DagMC* DAG;
}

//==============================================================================
// Non-member functions
//==============================================================================

void load_dagmc_geometry();
void free_memory_dagmc();
void read_geometry_dagmc();
bool read_uwuw_materials(pugi::xml_document& doc);
bool get_uwuw_materials_xml(std::string& s);

} // namespace openmc

#endif // DAGMC

#endif // OPENMC_DAGMC_H
