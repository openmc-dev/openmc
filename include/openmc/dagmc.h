
#ifndef OPENMC_DAGMC_H
#define OPENMC_DAGMC_H

namespace openmc {
extern "C" const bool DAGMC_ENABLED;
}

#ifdef DAGMC

#include "DagMC.hpp"

#include "openmc/position.h"
#include "openmc/xml_interface.h"

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
void read_dagmc_universes(pugi::xml_node node);
void read_dagmc_materials();
bool read_uwuw_materials(pugi::xml_document& doc);
bool get_uwuw_materials_xml(std::string& s);

} // namespace openmc

#else

inline void read_dagmc_materials() {};

#endif

#endif // OPENMC_DAGMC_H
