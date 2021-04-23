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

  //! Initialize the DAGMC accel. data structures, indices, material assignments, etc.
  void initialize();

  //! Reads UWUW materials and returns an ID map
  void read_uwuw_materials();
  //! Indicates whether or not UWUW materials are present
  //! \return True if UWUW materials are present, False if not
  bool uses_uwuw() const;

  //! Returns the index to the implicit complement's index in OpenMC for this DAGMC universe
  int32_t implicit_complement_idx() const;

  //! Transform UWUW materials into an OpenMC-readable XML format
  //! \return A string representing a materials.xml file of the UWUW materials in this universe
  std::string get_uwuw_materials_xml() const;

  //! Writes the UWUW material file to XML (for debugging purposes)
  void write_uwuw_materials_xml(const std::string& outfile = "uwuw_materials.xml") const;

  //! Assign a material to a cell based
  //! \param[in] mat_string The DAGMC material assignment string
  //! \param[in] c The OpenMC cell to which the material is assigned
  void legacy_assign_material(std::string mat_string,
                              std::unique_ptr<DAGCell>& c) const;

  //! Generate a string representing the ranges of IDs present in the DAGMC model
  //! \param[in] dim Dimension of the entities
  //! \return A string of the ID ranges for entities of dimension \p dim
  std::string dagmc_ids_for_dim(int dim) const;

  virtual bool find_cell(Particle &p) const override;

  // Data Members
  std::string filename_; //!< Name of the DAGMC file used to create this universe
  std::shared_ptr<moab::DagMC> dagmc_instance_; //!< DAGMC Instance for this universe
  std::shared_ptr<UWUW> uwuw_; //!< Pointer to the UWUW instance for this universe
  int32_t cell_idx_offset_; //!< An offset to the start of the cells in this universe in OpenMC's cell vector
  int32_t surf_idx_offset_; //!< An offset to the start of the surfaces in this universe in OpenMC's surface vector
  bool adjust_geometry_ids_; //!< Indicates whether or not to automatically generate new cell and surface IDs for the universe
  bool adjust_material_ids_; //!< Indicates whether or not to automatically generate new material IDs for the universe
};

//==============================================================================
// Non-member functions
//==============================================================================

void read_dagmc_universes(pugi::xml_node node);
int32_t next_cell(DAGUniverse* dag_univ, DAGCell* cur_cell, DAGSurface* surf_xed);

} // namespace openmc

#endif // DAGMC

#endif // OPENMC_DAGMC_H
