#ifndef OPENMC_DAGMC_H
#define OPENMC_DAGMC_H

namespace openmc {
extern "C" const bool DAGMC_ENABLED;
extern "C" const bool UWUW_ENABLED;
} // namespace openmc

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
#include "dagmcmetadata.hpp"

#include "openmc/cell.h"
#include "openmc/particle.h"
#include "openmc/position.h"
#include "openmc/surface.h"

class UWUW;

namespace openmc {

class DAGSurface : public Surface {
public:
  DAGSurface(std::shared_ptr<moab::DagMC> dag_ptr, int32_t dag_idx);

  moab::EntityHandle mesh_handle() const;

  double evaluate(Position r) const override;
  double distance(Position r, Direction u, bool coincident) const override;
  Direction normal(Position r) const override;
  Direction reflect(Position r, Direction u, GeometryState* p) const override;

  inline void to_hdf5_inner(hid_t group_id) const override {};

  // Accessor methods
  moab::DagMC* dagmc_ptr() const { return dagmc_ptr_.get(); }
  int32_t dag_index() const { return dag_index_; }

private:
  std::shared_ptr<moab::DagMC> dagmc_ptr_; //!< Pointer to DagMC instance
  int32_t dag_index_;                      //!< DagMC index of surface
};

class DAGCell : public Cell {
public:
  DAGCell(std::shared_ptr<moab::DagMC> dag_ptr, int32_t dag_idx);

  moab::EntityHandle mesh_handle() const;

  bool contains(Position r, Direction u, int32_t on_surface) const override;

  std::pair<double, int32_t> distance(Position r, Direction u,
    int32_t on_surface, GeometryState* p) const override;

  BoundingBox bounding_box() const override;

  void to_hdf5_inner(hid_t group_id) const override;

  // Accessor methods
  moab::DagMC* dagmc_ptr() const { return dagmc_ptr_.get(); }
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

  //! Alternative DAGMC universe constructor for external DAGMC instance
  explicit DAGUniverse(std::shared_ptr<moab::DagMC> external_dagmc_ptr,
    const std::string& filename = "", bool auto_geom_ids = false,
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

  //! Assign a material to a cell from uwuw material library
  //! \param[in] vol_handle The DAGMC material assignment string
  //! \param[in] c The OpenMC cell to which the material is assigned
  void uwuw_assign_material(
    moab::EntityHandle vol_handle, std::unique_ptr<DAGCell>& c) const;

  //! Assign a material to a cell based
  //! \param[in] mat_string The DAGMC material assignment string
  //! \param[in] c The OpenMC cell to which the material is assigned
  void legacy_assign_material(
    std::string mat_string, std::unique_ptr<DAGCell>& c) const;

  //! Assign a material overriding normal assignement to a cell
  //! \param[in] c The OpenMC cell to which the material is assigned
  void override_assign_material(std::unique_ptr<DAGCell>& c) const;

  //! Return the index into the model cells vector for a given DAGMC volume
  //! handle in the universe
  //! \param[in] vol MOAB handle to the DAGMC volume set
  int32_t cell_index(moab::EntityHandle vol) const;

  //! Return the index into the model surfaces vector for a given DAGMC surface
  //! handle in the universe
  //! \param[in] surf MOAB handle to the DAGMC surface set
  int32_t surface_index(moab::EntityHandle surf) const;

  //! Generate a string representing the ranges of IDs present in the DAGMC
  //! model. Contiguous chunks of IDs are represented as a range (i.e. 1-10). If
  //! there is a single ID a chunk, it will be represented as a single number
  //! (i.e. 2, 4, 6, 8). \param[in] dim Dimension of the entities \return A
  //! string of the ID ranges for entities of dimension \p dim
  std::string dagmc_ids_for_dim(int dim) const;

  bool find_cell(GeometryState& p) const override;

  void to_hdf5(hid_t universes_group) const override;

  // Data Members
  std::shared_ptr<moab::DagMC>
    dagmc_instance_;        //!< DAGMC Instance for this universe
  int32_t cell_idx_offset_; //!< An offset to the start of the cells in this
                            //!< universe in OpenMC's cell vector
  int32_t surf_idx_offset_; //!< An offset to the start of the surfaces in this
                            //!< universe in OpenMC's surface vector

  // Accessors
  moab::DagMC* dagmc_ptr() const { return dagmc_instance_.get(); }
  bool has_graveyard() const { return has_graveyard_; }

private:
  void set_id();        //!< Deduce the universe id from model::universes
  void init_dagmc();    //!< Create and initialise DAGMC pointer
  void init_metadata(); //!< Create and initialise dagmcMetaData pointer
  void init_geometry(); //!< Create cells and surfaces from DAGMC entities

  std::string
    filename_; //!< Name of the DAGMC file used to create this universe
  std::shared_ptr<UWUW>
    uwuw_; //!< Pointer to the UWUW instance for this universe
  std::unique_ptr<dagmcMetaData> dmd_ptr; //! Pointer to DAGMC metadata object
  bool adjust_geometry_ids_; //!< Indicates whether or not to automatically
                             //!< generate new cell and surface IDs for the
                             //!< universe
  bool adjust_material_ids_; //!< Indicates whether or not to automatically
                             //!< generate new material IDs for the universe
  bool has_graveyard_; //!< Indicates if the DAGMC geometry has a "graveyard"
                       //!< volume
  std::map<int, vector<std::string>>
    material_overrides; ///!< Map of material overrides
                        ///!< keys correspond to the DAGMCCell id
                        ///!< values are a list of materials used
                        ///!< git for the override
};

//==============================================================================
// Non-member functions
//==============================================================================

int32_t next_cell(int32_t surf, int32_t curr_cell, int32_t univ);

} // namespace openmc

#endif // DAGMC

#endif // OPENMC_DAGMC_H
