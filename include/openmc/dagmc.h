#ifndef OPENMC_DAGMC_H
#define OPENMC_DAGMC_H

namespace openmc {
extern "C" const bool DAGMC_ENABLED;
}

// always include the XML interface header
#include "openmc/xml_interface.h"

//==============================================================================
// Functions that are always defined
//==============================================================================

namespace openmc {

void read_dagmc_universes(pugi::xml_node node);
void check_dagmc_root_univ();

} // namespace openmc

#ifdef DAGMC

#include "DagMC.hpp"

#include "openmc/cell.h"
#include "openmc/particle.h"
#include "openmc/position.h"
#include "openmc/surface.h"

class UWUW;

namespace openmc {

class DAGSurface : public Surface {
public:
  DAGSurface(std::shared_ptr<moab::DagMC> dag_ptr, int32_t dag_idx);

  double evaluate(Position r) const override;
  double distance(Position r, Direction u, bool coincident) const override;
  Direction normal(Position r) const override;
  Direction reflect(Position r, Direction u, Particle* p) const override;

  inline void to_hdf5_inner(hid_t group_id) const override {};

  // Accessor methods
  const std::shared_ptr<moab::DagMC>& dagmc_ptr() const { return dagmc_ptr_; }
  int32_t dag_index() const { return dag_index_; }

private:
  std::shared_ptr<moab::DagMC> dagmc_ptr_; //!< Pointer to DagMC instance
  int32_t dag_index_;                      //!< DagMC index of surface
};

class DAGCell : public Cell {
public:
  DAGCell(std::shared_ptr<moab::DagMC> dag_ptr, int32_t dag_idx);

  bool contains(Position r, Direction u, int32_t on_surface) const override;

  std::pair<double, int32_t> distance(
    Position r, Direction u, int32_t on_surface, Particle* p) const override;

  BoundingBox bounding_box() const override;

  void to_hdf5_inner(hid_t group_id) const override;

  // Accessor methods
  const std::shared_ptr<moab::DagMC>& dagmc_ptr() const { return dagmc_ptr_; }
  int32_t dag_index() const { return dag_index_; }

private:
  std::shared_ptr<moab::DagMC> dagmc_ptr_; //!< Pointer to DagMC instance
  int32_t dag_index_;                      //!< DagMC index of cell
};

class DAGUniverse : public Universe {

public:
  explicit DAGUniverse(pugi::xml_node node);

  //! Create a new DAGMC universe
  //! \param[in] filename Name of the DAGMC file
  //! \param[in] auto_geom_ids Whether or not to automatically assign cell and
  //! surface IDs \param[in] auto_mat_ids Whether or not to automatically assign
  //! material IDs
  explicit DAGUniverse(const std::string& filename, bool auto_geom_ids = false,
    bool auto_mat_ids = false);

  //! Initialize the DAGMC accel. data structures, indices, material
  //! assignments, etc.
  void initialize();

  //! Reads UWUW materials and returns an ID map
  void read_uwuw_materials();
  //! Indicates whether or not UWUW materials are present
  //! \return True if UWUW materials are present, False if not
  bool uses_uwuw() const;

  //! Returns the index to the implicit complement's index in OpenMC for this
  //! DAGMC universe
  int32_t implicit_complement_idx() const;

  //! Transform UWUW materials into an OpenMC-readable XML format
  //! \return A string representing a materials.xml file of the UWUW materials
  //! in this universe
  std::string get_uwuw_materials_xml() const;

  //! Writes the UWUW material file to XML (for debugging purposes)
  void write_uwuw_materials_xml(
    const std::string& outfile = "uwuw_materials.xml") const;

  //! Assign a material to a cell based
  //! \param[in] mat_string The DAGMC material assignment string
  //! \param[in] c The OpenMC cell to which the material is assigned
  void legacy_assign_material(
    std::string mat_string, std::unique_ptr<DAGCell>& c) const;

  //! Generate a string representing the ranges of IDs present in the DAGMC
  //! model. Contiguous chunks of IDs are represented as a range (i.e. 1-10). If
  //! there is a single ID a chunk, it will be represented as a single number
  //! (i.e. 2, 4, 6, 8). \param[in] dim Dimension of the entities \return A
  //! string of the ID ranges for entities of dimension \p dim
  std::string dagmc_ids_for_dim(int dim) const;

  bool find_cell(Particle& p) const override;

  void to_hdf5(hid_t universes_group) const override;

  // Data Members
  std::shared_ptr<moab::DagMC>
    dagmc_instance_;        //!< DAGMC Instance for this universe
  int32_t cell_idx_offset_; //!< An offset to the start of the cells in this
                            //!< universe in OpenMC's cell vector
  int32_t surf_idx_offset_; //!< An offset to the start of the surfaces in this
                            //!< universe in OpenMC's surface vector

  // Accessors
  bool has_graveyard() const { return has_graveyard_; }

private:
  std::string
    filename_; //!< Name of the DAGMC file used to create this universe
  std::shared_ptr<UWUW>
    uwuw_;                   //!< Pointer to the UWUW instance for this universe
  bool adjust_geometry_ids_; //!< Indicates whether or not to automatically
                             //!< generate new cell and surface IDs for the
                             //!< universe
  bool adjust_material_ids_; //!< Indicates whether or not to automatically
                             //!< generate new material IDs for the universe
  bool has_graveyard_; //!< Indicates if the DAGMC geometry has a "graveyard"
                       //!< volume
};

//==============================================================================
// Non-member functions
//==============================================================================

int32_t next_cell(
  DAGUniverse* dag_univ, DAGCell* cur_cell, DAGSurface* surf_xed);

} // namespace openmc

#endif // DAGMC

#endif // OPENMC_DAGMC_H
