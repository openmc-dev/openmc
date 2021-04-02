#ifndef OPENMC_DAGMC_H
#define OPENMC_DAGMC_H

namespace openmc {
extern "C" const bool DAGMC_ENABLED;
}

#ifdef DAGMC

#include "DagMC.hpp"

#include "openmc/cell.h"
#include "openmc/particle.h"
#include "openmc/position.h"
#include "openmc/surface.h"
#include "openmc/xml_interface.h"

class UWUW;

namespace openmc {

class DAGSurface : public Surface
{
public:
  DAGSurface();

  double evaluate(Position r) const;
  double distance(Position r, Direction u, bool coincident) const;
  Direction normal(Position r) const;
  Direction reflect(Position r, Direction u, Particle* p) const;

  inline void to_hdf5_inner(hid_t group_id) const override {};

  std::shared_ptr<moab::DagMC> dagmc_ptr_; //!< Pointer to DagMC instance
  int32_t dag_index_;      //!< DagMC index of surface
};

class DAGCell : public Cell {
public:
  DAGCell();

  bool contains(Position r, Direction u, int32_t on_surface) const;

  std::pair<double, int32_t>
  distance(Position r, Direction u, int32_t on_surface, Particle* p) const;

  BoundingBox bounding_box() const;

  void to_hdf5_inner(hid_t group_id) const;

  std::shared_ptr<moab::DagMC> dagmc_ptr_; //!< Pointer to DagMC instance
  int32_t dag_index_;      //!< DagMC index of cell
};

class DAGUniverse : public Universe {

public:
  explicit DAGUniverse(pugi::xml_node node);
  explicit DAGUniverse(const std::string& filename, bool auto_geom_ids = false);

  void initialize(); //!< Sets up the DAGMC instance and OpenMC internals

  void read_uwuw_materials(); //!< Reads UWUW materials and returns an ID map
  bool uses_uwuw() const;

  int32_t implicit_complement_idx() const;

  std::string get_uwuw_materials_xml() const;

  void write_uwuw_materials_xml(const std::string& outfile = "uwuw_materials.xml") const;

  void legacy_assign_material(std::string mat_string,
                              std::unique_ptr<DAGCell>& c) const;

  std::string dagmc_ids_for_dim(int dim) const;

  virtual bool find_cell(Particle &p) const override;

  // Data Members
  std::string filename_;
  std::shared_ptr<moab::DagMC> dagmc_instance_; //! DAGMC Instance for this universe
  std::shared_ptr<UWUW> uwuw_;
  int32_t cell_idx_offset_;
  int32_t surf_idx_offset_;
  bool adjust_geometry_ids_;
  bool adjust_material_ids_;
};

// DAGMC Functions
void read_dagmc_universes(pugi::xml_node node);
int32_t next_cell(DAGUniverse* dag_univ, DAGCell* cur_cell, DAGSurface* surf_xed);

} // namespace openmc

#endif // DAGMC

#endif // OPENMC_DAGMC_H
