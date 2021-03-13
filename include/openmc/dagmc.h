
#ifndef OPENMC_DAGMC_H
#define OPENMC_DAGMC_H

namespace openmc {
extern "C" const bool DAGMC_ENABLED;
}

#ifdef DAGMC

#include "DagMC.hpp"

#include "openmc/cell.h"
#include "openmc/position.h"
#include "openmc/xml_interface.h"

class UWUW;

namespace openmc {

namespace model {
  extern std::shared_ptr<moab::DagMC> DAG;
}


class DAGUniverse : public Universe {

public:
  explicit DAGUniverse(pugi::xml_node node);
  explicit DAGUniverse(const std::string& filename, bool auto_geom_ids = false);

  void initialize(); //!< Sets up the DAGMC instance and OpenMC internals

  std::shared_ptr<UWUW>
  read_uwuw_materials(); //!< Reads UWUW materials and returns an ID map

  // Data Members
  std::string filename_;
  std::shared_ptr<moab::DagMC> dagmc_instance_; //! DAGMC Instance for this universe
  int32_t cell_idx_offset_;
  int32_t surf_idx_offset_;
  bool adjust_geometry_ids_;
  bool adjust_material_ids_;
};

void load_dagmc_geometry();
int32_t create_dagmc_universe(const std::string& filename);
void read_geometry_dagmc();
void read_dagmc_universes(pugi::xml_node node);
bool read_uwuw_materials(pugi::xml_document& doc);
bool get_uwuw_materials_xml(std::string& s);

} // namespace openmc

#else

inline void read_dagmc_materials() {};

#endif

#endif // OPENMC_DAGMC_H
