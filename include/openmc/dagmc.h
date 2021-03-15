
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

class DAGUniverse : public Universe {

public:
  explicit DAGUniverse(pugi::xml_node node);
  explicit DAGUniverse(const std::string& filename, bool auto_geom_ids = false);

  void initialize(); //!< Sets up the DAGMC instance and OpenMC internals

  void read_uwuw_materials(); //!< Reads UWUW materials and returns an ID map

  bool uses_uwuw() const;

  std::string get_uwuw_materials_xml() const;

  void write_uwuw_materials_xml(const std::string& outfile = "uwuw_materials.xml") const;

  void legacy_assign_material(std::string mat_string,
                              std::unique_ptr<DAGCell>& c) const;

  std::string dagmc_ids_for_dim(int dim) const;

  // Data Members
  std::string filename_;
  std::shared_ptr<moab::DagMC> dagmc_instance_; //! DAGMC Instance for this universe
  std::shared_ptr<UWUW> uwuw_;
  int32_t cell_idx_offset_;
  int32_t surf_idx_offset_;
  bool adjust_geometry_ids_;
  bool adjust_material_ids_;
};

void read_dagmc_universes(pugi::xml_node node);

} // namespace openmc

#else

inline void read_dagmc_materials() {};

#endif

#endif // OPENMC_DAGMC_H
