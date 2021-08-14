//! \file mesh.h
//! \brief Mesh types used for tallies, Shannon entropy, CMFD, etc.

#ifndef OPENMC_MESH_H
#define OPENMC_MESH_H

#include <unordered_map>

#include "hdf5.h"
#include "pugixml.hpp"
#include "xtensor/xtensor.hpp"

#include "openmc/memory.h" // for unique_ptr
#include "openmc/particle.h"
#include "openmc/position.h"
#include "openmc/vector.h"

#ifdef DAGMC
#include "moab/AdaptiveKDTree.hpp"
#include "moab/Core.hpp"
#include "moab/GeomUtil.hpp"
#include "moab/Matrix3.hpp"
#endif

#ifdef LIBMESH
#include "libmesh/bounding_box.h"
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/equation_systems.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/explicit_system.h"
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/point.h"
#endif

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

extern "C" const bool LIBMESH_ENABLED;

class Mesh;

namespace model {

extern std::unordered_map<int32_t, int32_t> mesh_map;
extern vector<unique_ptr<Mesh>> meshes;

} // namespace model

#ifdef LIBMESH
namespace settings {
// used when creating new libMesh::Mesh instances
extern unique_ptr<libMesh::LibMeshInit> libmesh_init;
extern const libMesh::Parallel::Communicator* libmesh_comm;
} // namespace settings
#endif

class Mesh {
public:
  // Constructors and destructor
  Mesh() = default;
  Mesh(pugi::xml_node node);
  virtual ~Mesh() = default;

  // Methods

  //! Determine which bins were crossed by a particle
  //
  //! \param[in] r0 Previous position of the particle
  //! \param[in] r1 Current position of the particle
  //! \param[in] u Particle direction
  //! \param[out] bins Bins that were crossed
  //! \param[out] lengths Fraction of tracklength in each bin
  virtual void bins_crossed(Position r0, Position r1, const Direction& u,
    vector<int>& bins, vector<double>& lengths) const = 0;

  //! Determine which surface bins were crossed by a particle
  //
  //! \param[in] r0 Previous position of the particle
  //! \param[in] r1 Current position of the particle
  //! \param[in] u Particle direction
  //! \param[out] bins Surface bins that were crossed
  virtual void surface_bins_crossed(
    Position r0, Position r1, const Direction& u, vector<int>& bins) const = 0;

  //! Get bin at a given position in space
  //
  //! \param[in] r Position to get bin for
  //! \return Mesh bin
  virtual int get_bin(Position r) const = 0;

  //! Get the number of mesh cells.
  virtual int n_bins() const = 0;

  //! Get the number of mesh cell surfaces.
  virtual int n_surface_bins() const = 0;

  //! Set the mesh ID
  void set_id(int32_t id = -1);

  //! Write mesh data to an HDF5 group
  //
  //! \param[in] group HDF5 group
  virtual void to_hdf5(hid_t group) const = 0;

  //! Find the mesh lines that intersect an axis-aligned slice plot
  //
  //! \param[in] plot_ll The lower-left coordinates of the slice plot.
  //! \param[in] plot_ur The upper-right coordinates of the slice plot.
  //! \return A pair of vectors indicating where the mesh lines lie along each
  //!   of the plot's axes.  For example an xy-slice plot will get back a vector
  //!   of x-coordinates and another of y-coordinates.  These vectors may be
  //!   empty for low-dimensional meshes.
  virtual std::pair<vector<double>, vector<double>> plot(
    Position plot_ll, Position plot_ur) const = 0;

  //! Return a string representation of the mesh bin
  //
  //! \param[in] bin Mesh bin to generate a label for
  virtual std::string bin_label(int bin) const = 0;

  // Data members
  int id_ {-1};     //!< User-specified ID
  int n_dimension_; //!< Number of dimensions
};

class StructuredMesh : public Mesh {
public:
  StructuredMesh() = default;
  StructuredMesh(pugi::xml_node node) : Mesh {node} {};
  virtual ~StructuredMesh() = default;

  int get_bin(Position r) const override;

  int n_bins() const override;

  int n_surface_bins() const override;

  void bins_crossed(Position r0, Position r1, const Direction& u,
    vector<int>& bins, vector<double>& lengths) const override;

  //! Count number of bank sites in each mesh bin / energy bin
  //
  //! \param[in] Pointer to bank sites
  //! \param[in] Number of bank sites
  //! \param[out] Whether any bank sites are outside the mesh
  xt::xtensor<double, 1> count_sites(
    const SourceSite* bank, int64_t length, bool* outside) const;

  //! Get bin given mesh indices
  //
  //! \param[in] Array of mesh indices
  //! \return Mesh bin
  virtual int get_bin_from_indices(const int* ijk) const;

  //! Get mesh indices given a position
  //
  //! \param[in] r Position to get indices for
  //! \param[out] ijk Array of mesh indices
  //! \param[out] in_mesh Whether position is in mesh
  virtual void get_indices(Position r, int* ijk, bool* in_mesh) const;

  //! Get mesh indices corresponding to a mesh bin
  //
  //! \param[in] bin Mesh bin
  //! \param[out] ijk Mesh indices
  virtual void get_indices_from_bin(int bin, int* ijk) const;

  //! Get mesh index in a particular direction
  //!
  //! \param[in] r Coordinate to get index for
  //! \param[in] i Direction index
  virtual int get_index_in_direction(double r, int i) const = 0;

  //! Check where a line segment intersects the mesh and if it intersects at all
  //
  //! \param[in,out] r0 In: starting position, out: intersection point
  //! \param[in] r1 Ending position
  //! \param[out] ijk Indices of the mesh bin containing the intersection point
  //! \return Whether the line segment connecting r0 and r1 intersects mesh
  virtual bool intersects(Position& r0, Position r1, int* ijk) const;

  //! Get the coordinate for the mesh grid boundary in the positive direction
  //!
  //! \param[in] ijk Array of mesh indices
  //! \param[in] i Direction index
  virtual double positive_grid_boundary(int* ijk, int i) const = 0;

  //! Get the coordinate for the mesh grid boundary in the negative direction
  //!
  //! \param[in] ijk Array of mesh indices
  //! \param[in] i Direction index
  virtual double negative_grid_boundary(int* ijk, int i) const = 0;

  //! Get a label for the mesh bin
  std::string bin_label(int bin) const override;

  // Data members
  xt::xtensor<double, 1> lower_left_;  //!< Lower-left coordinates of mesh
  xt::xtensor<double, 1> upper_right_; //!< Upper-right coordinates of mesh
  xt::xtensor<int, 1> shape_; //!< Number of mesh elements in each dimension

protected:
  virtual bool intersects_1d(Position& r0, Position r1, int* ijk) const;
  virtual bool intersects_2d(Position& r0, Position r1, int* ijk) const;
  virtual bool intersects_3d(Position& r0, Position r1, int* ijk) const;
};

//==============================================================================
//! Tessellation of n-dimensional Euclidean space by congruent squares or cubes
//==============================================================================

class RegularMesh : public StructuredMesh {
public:
  // Constructors
  RegularMesh() = default;
  RegularMesh(pugi::xml_node node);

  // Overridden methods
  void surface_bins_crossed(Position r0, Position r1, const Direction& u,
    vector<int>& bins) const override;

  int get_index_in_direction(double r, int i) const override;

  double positive_grid_boundary(int* ijk, int i) const override;

  double negative_grid_boundary(int* ijk, int i) const override;

  std::pair<vector<double>, vector<double>> plot(
    Position plot_ll, Position plot_ur) const override;

  void to_hdf5(hid_t group) const override;

  // New methods

  //! Count number of bank sites in each mesh bin / energy bin
  //
  //! \param[in] bank Array of bank sites
  //! \param[out] Whether any bank sites are outside the mesh
  //! \return Array indicating number of sites in each mesh/energy bin
  xt::xtensor<double, 1> count_sites(
    const SourceSite* bank, int64_t length, bool* outside) const;

  // Data members
  double volume_frac_;           //!< Volume fraction of each mesh element
  xt::xtensor<double, 1> width_; //!< Width of each mesh element
};

class RectilinearMesh : public StructuredMesh {
public:
  // Constructors
  RectilinearMesh() = default;
  RectilinearMesh(pugi::xml_node node);

  // Overridden methods
  void surface_bins_crossed(Position r0, Position r1, const Direction& u,
    vector<int>& bins) const override;

  int get_index_in_direction(double r, int i) const override;

  double positive_grid_boundary(int* ijk, int i) const override;

  double negative_grid_boundary(int* ijk, int i) const override;

  std::pair<vector<double>, vector<double>> plot(
    Position plot_ll, Position plot_ur) const override;

  void to_hdf5(hid_t group) const override;

  vector<vector<double>> grid_;

  int set_grid();
};

// Abstract class for unstructured meshes
class UnstructuredMesh : public Mesh {

public:
  // Constructors
  UnstructuredMesh() {};
  UnstructuredMesh(pugi::xml_node node);
  UnstructuredMesh(const std::string& filename);

  // Overridden Methods
  void surface_bins_crossed(Position r0, Position r1, const Direction& u,
    vector<int>& bins) const override;

  void to_hdf5(hid_t group) const override;

  std::string bin_label(int bin) const override;

  // Methods

  //! Add a variable to the mesh instance
  virtual void add_score(const std::string& var_name) = 0;

  //! Remove tally data from the instance
  virtual void remove_scores() = 0;

  //! Set the value of a bin for a variable on the internal
  //  mesh instance
  virtual void set_score_data(const std::string& var_name,
    const vector<double>& values, const vector<double>& std_dev) = 0;

  //! Write the unstructured mesh to file
  //
  //! \param[in] filename Base of the file to write
  virtual void write(const std::string& base_filename) const = 0;

  //! Retrieve a centroid for the mesh cell
  //
  //! \param[in] bin Bin to return the centroid for
  //! \return The centroid of the bin
  virtual Position centroid(int bin) const = 0;

  //! Get the volume of a mesh bin
  //
  //! \param[in] bin Bin to return the volume for
  //! \return Volume of the bin
  virtual double volume(int bin) const = 0;

  //! Get the library used for this unstructured mesh
  virtual std::string library() const = 0;

  // Data members
  bool output_ {
    true}; //!< Write tallies onto the unstructured mesh at the end of a run
  std::string filename_; //!< Path to unstructured mesh file

private:
  //! Setup method for the mesh. Builds data structures,
  //! sets up element mapping, creates bounding boxes, etc.
  virtual void initialize() = 0;
};

#ifdef DAGMC

class MOABMesh : public UnstructuredMesh {
public:
  // Constructors
  MOABMesh() = default;
  MOABMesh(pugi::xml_node);
  MOABMesh(const std::string& filename);
  MOABMesh(std::shared_ptr<moab::Interface> external_mbi);

  // Overridden Methods

  void bins_crossed(Position r0, Position r1, const Direction& u,
    vector<int>& bins, vector<double>& lengths) const override;

  int get_bin(Position r) const override;

  int n_bins() const override;

  int n_surface_bins() const override;

  std::pair<vector<double>, vector<double>> plot(
    Position plot_ll, Position plot_ur) const override;

  std::string library() const override;

  //! Add a score to the mesh instance
  void add_score(const std::string& score) override;

  //! Remove all scores from the mesh instance
  void remove_scores() override;

  //! Set data for a score
  void set_score_data(const std::string& score, const vector<double>& values,
    const vector<double>& std_dev) override;

  //! Write the mesh with any current tally data
  void write(const std::string& base_filename) const override;

  Position centroid(int bin) const override;

  double volume(int bin) const override;

private:
  void initialize() override;

  // Methods

  //! Create the MOAB interface pointer
  void create_interface();

  //! Find all intersections with faces of the mesh.
  //
  //! \param[in] start Staring location
  //! \param[in] dir Normalized particle direction
  //! \param[in] track_len length of particle track
  //! \param[out] Mesh intersections
  void intersect_track(const moab::CartVect& start, const moab::CartVect& dir,
    double track_len, vector<double>& hits) const;

  //! Calculate the volume for a given tetrahedron handle.
  //
  // \param[in] tet MOAB EntityHandle of the tetrahedron
  double tet_volume(moab::EntityHandle tet) const;

  //! Find the tetrahedron for the given location if
  //! one exists
  //
  //! \param[in]
  //! \return MOAB EntityHandle of tet
  moab::EntityHandle get_tet(const Position& r) const;

  //! Return the containing tet given a position
  moab::EntityHandle get_tet(const moab::CartVect& r) const
  {
    return get_tet(Position(r[0], r[1], r[2]));
  };

  //! Check for point containment within a tet; uses
  //! pre-computed barycentric data.
  //
  //! \param[in] r Position to check
  //! \param[in] MOAB terahedron to check
  //! \return True if r is inside, False if r is outside
  bool point_in_tet(const moab::CartVect& r, moab::EntityHandle tet) const;

  //! Compute barycentric coordinate data for all tetrahedra
  //! in the mesh.
  //
  //! \param[in] tets MOAB Range of tetrahedral elements
  void compute_barycentric_data(const moab::Range& tets);

  //! Translate a MOAB EntityHandle to its corresponding bin.
  //
  //! \param[in] eh MOAB EntityHandle to translate
  //! \return Mesh bin
  int get_bin_from_ent_handle(moab::EntityHandle eh) const;

  //! Translate a bin to its corresponding MOAB EntityHandle
  //! for the tetrahedron representing that bin.
  //
  //! \param[in] bin Bin value to translate
  //! \return MOAB EntityHandle of tet
  moab::EntityHandle get_ent_handle_from_bin(int bin) const;

  //! Get the bin for a given mesh cell index
  //
  //! \param[in] idx Index of the mesh cell.
  //! \return Mesh bin
  int get_bin_from_index(int idx) const;

  //! Get the mesh cell index for a given position
  //
  //! \param[in] r Position to get index for
  //! \param[in,out] in_mesh Whether position is in the mesh
  int get_index(const Position& r, bool* in_mesh) const;

  //! Get the mesh cell index from a bin
  //
  //! \param[in] bin Bin to get the index for
  //! \return Index of the bin
  int get_index_from_bin(int bin) const;

  //! Build a KDTree for all tetrahedra in the mesh. All
  //! triangles representing 2D faces of the mesh are
  //! added to the tree as well.
  //
  //! \param[in] all_tets MOAB Range of tetrahedra for the tree
  void build_kdtree(const moab::Range& all_tets);

  //! Get the tags for a score from the mesh instance
  //! or create them if they are not there
  //
  //! \param[in] score Name of the score
  //! \return The MOAB value and error tag handles, respectively
  std::pair<moab::Tag, moab::Tag> get_score_tags(std::string score) const;

  // Data members
  moab::Range ehs_; //!< Range of tetrahedra EntityHandle's in the mesh
  moab::EntityHandle tetset_;      //!< EntitySet containing all tetrahedra
  moab::EntityHandle kdtree_root_; //!< Root of the MOAB KDTree
  std::shared_ptr<moab::Interface> mbi_;    //!< MOAB instance
  unique_ptr<moab::AdaptiveKDTree> kdtree_; //!< MOAB KDTree instance
  vector<moab::Matrix3> baryc_data_;        //!< Barycentric data for tetrahedra
  vector<std::string> tag_names_; //!< Names of score tags added to the mesh
};

#endif

#ifdef LIBMESH

class LibMesh : public UnstructuredMesh {
public:
  // Constructors
  LibMesh(pugi::xml_node node);
  LibMesh(const std::string& filename);

  // Overridden Methods
  void bins_crossed(Position r0, Position r1, const Direction& u,
    vector<int>& bins, vector<double>& lengths) const override;

  int get_bin(Position r) const override;

  int n_bins() const override;

  int n_surface_bins() const override;

  std::pair<vector<double>, vector<double>> plot(
    Position plot_ll, Position plot_ur) const override;

  std::string library() const override;

  void add_score(const std::string& var_name) override;

  void remove_scores() override;

  void set_score_data(const std::string& var_name, const vector<double>& values,
    const vector<double>& std_dev) override;

  void write(const std::string& base_filename) const override;

  Position centroid(int bin) const override;

  double volume(int bin) const override;

private:
  void initialize() override;

  // Methods

  //! Translate a bin value to an element reference
  const libMesh::Elem& get_element_from_bin(int bin) const;

  //! Translate an element pointer to a bin index
  int get_bin_from_element(const libMesh::Elem* elem) const;

  // Data members
  unique_ptr<libMesh::Mesh> m_; //!< pointer to the libMesh mesh instance
  vector<unique_ptr<libMesh::PointLocatorBase>>
    pl_; //!< per-thread point locators
  unique_ptr<libMesh::EquationSystems>
    equation_systems_; //!< pointer to the equation systems of the mesh
  std::string
    eq_system_name_; //!< name of the equation system holding OpenMC results
  std::unordered_map<std::string, unsigned int>
    variable_map_; //!< mapping of variable names (tally scores) to libMesh
                   //!< variable numbers
  libMesh::BoundingBox bbox_; //!< bounding box of the mesh
  libMesh::dof_id_type
    first_element_id_; //!< id of the first element in the mesh
};

#endif

//==============================================================================
// Non-member functions
//==============================================================================

//! Read meshes from either settings/tallies
//
//! \param[in] root XML node
void read_meshes(pugi::xml_node root);

//! Write mesh data to an HDF5 group
//
//! \param[in] group HDF5 group
void meshes_to_hdf5(hid_t group);

void free_memory_mesh();

} // namespace openmc

#endif // OPENMC_MESH_H
