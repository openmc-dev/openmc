#ifndef OPENMC_XDG_H
#define OPENMC_XDG_H

namespace openmc {
extern "C" const bool XDG_ENABLED;
}

// always include the XML interface header
#include "openmc/xml_interface.h"

namespace openmc {

void read_xdg_universes(pugi::xml_node node);
void read_xdg_mesh_universes(pugi::xml_node node);

} // namespace openmc

#ifdef OPENMC_XDG_ENABLED

#include "xdg/xdg.h"

#include "openmc/cell.h"
#include "openmc/mesh.h"
#include "openmc/particle.h"
#include "openmc/position.h"
#include "openmc/surface.h"
namespace openmc {

class XDGMesh : public UnstructuredMesh {

public:
  //----------------------------------------------------------------------------
  // Constructors
  XDGMesh() = default;
  XDGMesh(pugi::xml_node node);
  XDGMesh(hid_t group);
  XDGMesh(const std::string& filename, double length_multiplier = 1.0);
  XDGMesh(
    std::shared_ptr<xdg::XDG> external_xdg, double length_multiplier = 1.0);

  static const std::string mesh_interface;

  //----------------------------------------------------------------------------
  // Methods

  //! Get the underlying XDG instance
  //!
  //! \return Shared pointer to the XDG instance
  const std::shared_ptr<xdg::XDG>& xdg_instance() const { return xdg_; }

  //! Check whether a bin index is valid
  //!
  //! \param[in] bin Bin index to check
  //! \return True if the bin index is in [0, n_bins())
  bool bin_is_valid(int bin) const { return bin >= 0 && bin < n_bins(); }

  //! Convert a mesh bin index to an XDG MeshID
  //!
  //! \param[in] bin Bin index to convert
  //! \return XDG MeshID corresponding to the bin
  xdg::MeshID bin_to_mesh_id(int bin) const;

  //! Convert an XDG MeshID to a mesh bin index
  //!
  //! \param[in] id XDG MeshID to convert
  //! \return Bin index corresponding to the MeshID, or -1 if invalid
  int mesh_id_to_bin(xdg::MeshID id) const;

  //! Get the name of the underlying mesh library
  //!
  //! \return Name of the mesh library being used in XDG
  std::string library() const override;

  //----------------------------------------------------------------------------
  // Overridden Methods

  //! Perform any preparation needed to support use in mesh filters
  void prepare_for_point_location() override;

  Position sample_element(int32_t bin, uint64_t* seed) const override;

  void bins_crossed(Position r0, Position r1, const Direction& u,
    vector<int>& bins, vector<double>& lengths) const override;

  int get_bin(Position r) const override;

  int n_bins() const override;

  int n_surface_bins() const override;

  std::pair<vector<double>, vector<double>> plot(
    Position plot_ll, Position plot_ur) const override;

  //! Add a score to the mesh instance
  void add_score(const std::string& score) override {};

  //! Remove all scores from the mesh instance
  void remove_scores() override {};

  //! Set data for a score
  void set_score_data(const std::string& score, const vector<double>& values,
    const vector<double>& std_dev) override {};

  //! Write the mesh with any current tally data
  void write(const std::string& base_filename) const override;

  Position centroid(int bin) const override;

  int n_vertices() const override;

  Position vertex(int id) const override;

  std::vector<int> connectivity(int id) const override;

  //! Get the volume of a mesh bin
  //
  //! \param[in] bin Bin to return the volume for
  //! \return Volume of the bin
  double volume(int bin) const override;

  //! Get the distance to the nearest boundary for a given position and
  //! direction \param[in] g GeometryState object containing position and
  //! direction \return NextMeshCell struct containing distance, face index, and
  //! next indices
  NextMeshCell distance_to_bin_boundary(GeometryState& g) const;

  //! Get the distance to the nearest boundary for a given position and
  //! direction \param[in] r Position to check \param[in] u Direction to check
  //! \return Distance to the nearest boundary
  NextMeshCell distance_to_bin_boundary(
    int bin, const Position& r, const Direction& u) const;

private:
  void initialize() override;

  int n_surface_bins_ {0}; //!< Number of surface bins in the mesh

  //----------------------------------------------------------------------------
  // Private data members
  std::shared_ptr<xdg::XDG> xdg_; //!< XDG instance
  xdg::MeshLibrary mesh_library_ {
    xdg::MeshLibrary::LIBMESH}; //!< Mesh library type
};

class XDGGeometryObject {
public:
  XDGGeometryObject(std::shared_ptr<xdg::XDG> xdg_ptr, xdg::MeshID xdg_id)
    : xdg_ptr_ {xdg_ptr}, xdg_id_ {xdg_id}
  {}

  // Accessor methods
  xdg::MeshID xdg_id() const { return xdg_id_; }
  const std::shared_ptr<xdg::XDG>& xdg_ptr() const { return xdg_ptr_; }

private:
  std::shared_ptr<xdg::XDG> xdg_ptr_; //!< Pointer to XDG instance
  xdg::MeshID xdg_id_;                //!< XDG ID
};

class XDGSurface : public XDGGeometryObject, public Surface {
public:
  XDGSurface(std::shared_ptr<xdg::XDG> dag_ptr, xdg::MeshID xdg_id);

  GeometryType geom_type() const override
  {
    return GeometryType::XDG_SURFACE_MESH;
  }

  double evaluate(Position r) const override;
  double distance(Position r, Direction u, bool coincident) const override;

  Direction normal(Position r) const override;
  Direction reflect(Position r, Direction u, GeometryState* p) const override;

  inline void to_hdf5_inner(hid_t group_id) const override {};
};

class XDGCell : public XDGGeometryObject, public Cell {
public:
  XDGCell(const std::shared_ptr<xdg::XDG>& xdg_ptr, xdg::MeshID xdg_id);

  GeometryType geom_type() const override
  {
    return GeometryType::XDG_SURFACE_MESH;
  }

  bool contains(Position r, Direction u, int32_t on_surface) const override;

  std::pair<double, int32_t> distance(Position r, Direction u,
    int32_t on_surface, GeometryState* p) const override;

  BoundingBox bounding_box() const override;

  void to_hdf5_inner(hid_t group_id) const override;
};

class XDGUniverse : public Universe {

public:
  explicit XDGUniverse(pugi::xml_node node);

  GeometryType geom_type() const override
  {
    return GeometryType::XDG_SURFACE_MESH;
  }

  //! Create a new XDG universe
  //! \param[in] filename Name of the XDG file
  //! \param[in] auto_geom_ids Whether or not to automatically assign cell and
  //! surface IDs
  //! \param[in] auto_mat_ids Whether or not to automatically assign
  //! material IDs
  explicit XDGUniverse(const std::string& filename, bool auto_geom_ids = false,
    bool auto_mat_ids = false);

  //! Alternative XDG universe constructor for external XDG instance
  explicit XDGUniverse(std::shared_ptr<xdg::XDG> external_xdg_ptr,
    const std::string& filename = "", bool auto_geom_ids = false,
    bool auto_mat_ids = false);

  //! Initialize the XDG accel. data structures, indices, material
  //! assignments, etc.
  void initialize();

  //! Returns the index to the implicit complement's index in OpenMC for this
  //! XDG universe
  int32_t implicit_complement_idx() const;

  //! Assign a material to a cell based
  //! \param[in] mat_string The XDG material assignment string
  //! \param[in] c The OpenMC cell to which the material is assigned
  void assign_material(
    std::string& mat_string, std::unique_ptr<XDGCell>& c) const;

  //! Return the index into the model cells vector for a given XDG volume
  //! handle in the universe
  //! \param[in] vol MOAB handle to the XDG volume set
  int32_t cell_index(xdg::MeshID voume) const;

  //! Return the index into the model surfaces vector for a given XDG surface
  //! handle in the universe
  //! \param[in] surf MOAB handle to the XDG surface set
  int32_t surface_index(xdg::MeshID surface) const;

  //! Generate a string representing the ranges of IDs present in the XDG
  //! model. Contiguous chunks of IDs are represented as a range (i.e. 1-10). If
  //! there is a single ID a chunk, it will be represented as a single number
  //! (i.e. 2, 4, 6, 8).
  //! \param[in] dim Dimension of the entities
  //! \return A string of the ID ranges for entities of dimension \p dim
  std::string xdg_ids_for_dim(int dim) const;

  //! Move the particle to the next cell in the XDG universe
  //! \param[in] p The particle to move. Particle should be on a
  //! surface when called.
  void next_cell(Particle& p) const;

  bool find_cell(GeometryState& p) const override;

  void to_hdf5(hid_t universes_group) const override;

  // Data Members
  std::shared_ptr<xdg::XDG> xdg_instance_; //!< XDG Instance for this universe
  int32_t cell_idx_offset_; //!< An offset to the start of the cells in this
                            //!< universe in OpenMC's cell vector
  int32_t surf_idx_offset_; //!< An offset to the start of the surfaces in this
                            //!< universe in OpenMC's surface vector

  std::string library() const { return xdg_mesh()->library(); }
  std::string filename() const { return xdg_mesh()->filename(); }

  // Accessors
  int32_t mesh_idx() const { return mesh_idx_; }
  const XDGMesh* xdg_mesh() const
  {
    return dynamic_cast<const XDGMesh*>(model::meshes[mesh_idx_].get());
  }
  const std::unique_ptr<Mesh>& mesh() const { return model::meshes[mesh_idx_]; }
  const std::shared_ptr<xdg::XDG>& xdg_ptr() const
  {
    return xdg_mesh()->xdg_instance();
  }

private:
  void set_id();        //!< Deduce the universe id from model::universes
  void init_xdg();      //!< Create and initialise XDG pointer
  void init_metadata(); //!< Create and initialise dagmcMetaData pointer
  void init_geometry(); //!< Create cells and surfaces from XDG entities

  bool adjust_geometry_ids_; //!< Indicates whether or not to automatically
                             //!< generate new cell and surface IDs for the
                             //!< universe
  bool adjust_material_ids_; //!< Indicates whether or not to automatically
                             //!< generate new material IDs for the universe

  int32_t mesh_idx_; //!< The index of the mesh in the model::meshes vector

  // mappings from XDG IDs to OpenMC surface and cell indices
  std::unordered_map<xdg::MeshID, int32_t> surface_index_map_;
  std::unordered_map<xdg::MeshID, int32_t> cell_index_map_;
};

class XDGMeshUniverse : public Universe {

public:
  // constructors
  XDGMeshUniverse() = default;

  GeometryType geom_type() const override
  {
    return GeometryType::XDG_VOLUME_MESH;
  }

  explicit XDGMeshUniverse(pugi::xml_node node);

  // setup functions

  // contains mesh-generic code
  void create_cells(pugi::xml_node node);

  // match a material name to an OpenMC material index
  int32_t match_material(const std::string& material_name) const;

  void set_boundary_conditions();

  // transport methods
  virtual bool find_cell(openmc::GeometryState& p) const override;

  void next_cell(Particle& p) const;

  // accessors
  int32_t outer_material() const { return outer_material_; }
  int32_t& outer_material() { return outer_material_; }

  const std::unordered_map<xdg::MeshID, std::vector<int32_t>>&
  element_material_map() const
  {
    return element_material_map_;
  }

  // Accessors
  int32_t mesh_idx() const { return mesh_idx_; }
  const std::unique_ptr<Mesh>& mesh() const { return model::meshes[mesh_idx_]; }
  const XDGMesh* xdg_mesh() const
  {
    return dynamic_cast<const XDGMesh*>(model::meshes[mesh_idx_].get());
  }
  const std::shared_ptr<xdg::XDG>& xdg_instance() const
  {
    return xdg_mesh()->xdg_instance();
  }

protected:
  int32_t mesh_idx_;
  std::string name_;
  std::unordered_map<xdg::MeshID, std::vector<int32_t>> element_material_map_;
  int32_t outer_material_ {MATERIAL_VOID};
};

class XDGMeshCell : public Cell {
public:
  XDGMeshCell(int32_t mesh, int32_t element_idx)
    : mesh_(mesh), elem_idx_(element_idx)
  {}

  GeometryType geom_type() const override
  {
    return GeometryType::XDG_VOLUME_MESH;
  }

  virtual bool contains(
    Position r, Direction u, int32_t on_surface) const override
  {
    int mesh_bin = model::meshes[mesh_]->get_bin(r);
    return mesh_bin == elem_idx_ - 1;
  };

  virtual std::pair<double, int32_t> distance(Position r, Direction u,
    int32_t on_surface, GeometryState* p) const override
  {
    // if this element is the background, determine
    // if the particle might re-enter the mesh
    if (elem_idx_ == C_NONE) {
      auto ipc = xdg_ptr()->mesh_manager()->implicit_complement();
      auto ipc_elem =
        xdg_ptr()->ray_fire(ipc, {r.x, r.y, r.z}, {u.x, u.y, u.z});
      if (ipc_elem.second == C_NONE) {
        return {INFTY, elem_idx_};
      }
      // if the particle will re-enter the mesh, return the distance to the
      // surface of the implicit complement
      auto new_r = r + u * (ipc_elem.first + TINY_BIT);
      auto next_element = xdg_ptr()->find_element({new_r.x, new_r.y, new_r.z});
      return {ipc_elem.first, next_element};
    }
    auto result = xdg_ptr()->mesh_manager()->next_element(
      elem_idx_, {r.x, r.y, r.z}, {u.x, u.y, u.z});
    return {result.second, result.first};
  }

  virtual int32_t material(int32_t instance) const override
  {
    return material_[0];
    // const auto& element_materials =
    // mesh_univ()->element_material_map().at(elem_idx_); return
    // element_materials.size() > 1 ? element_materials[instance] :
    // element_materials[0];
  }

  const XDGMeshUniverse* mesh_univ() const
  {
    return dynamic_cast<const XDGMeshUniverse*>(
      model::universes[universe_idx_].get());
  }

  const XDGMesh* xdg_mesh() const
  {
    return dynamic_cast<const XDGMesh*>(model::meshes[mesh_].get());
  }

  const xdg::XDG* xdg_ptr() const { return xdg_mesh()->xdg_instance().get(); }

  virtual void to_hdf5_inner(hid_t group_id) const override {};

  virtual BoundingBox bounding_box() const override { return BoundingBox {}; };

protected:
  int32_t universe_idx_;
  int32_t mesh_;
  int32_t elem_idx_;
};

//==============================================================================
// Non-member functions
//==============================================================================

int32_t xdg_next_cell(int32_t surf, int32_t curr_cell, int32_t univ);

} // namespace openmc

#endif // OPENMC_XDG_ENABLED

#endif // OPENMC_XDG_H
