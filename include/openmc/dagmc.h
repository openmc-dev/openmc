
#ifndef OPENMC_DAGMC_H
#define OPENMC_DAGMC_H

namespace openmc {
extern "C" const bool DAGMC_ENABLED;
}

#ifdef DAGMC

#include "DagMC.hpp"
#include "openmc/xml_interface.h"
#include "openmc/position.h"

namespace openmc {

namespace model {
  extern std::shared_ptr<moab::DagMC> DAG;
}

//==============================================================================
// Non-member functions
//==============================================================================

void load_dagmc_geometry();
int32_t create_dagmc_universe(const std::string& filename);
void read_geometry_dagmc();
bool read_uwuw_materials(pugi::xml_document& doc);
bool get_uwuw_materials_xml(std::string& s);

} // namespace openmc

#endif // DAGMC

#endif // OPENMC_DAGMC_H
