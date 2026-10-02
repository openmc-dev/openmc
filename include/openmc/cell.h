#ifndef OPENMC_CELL_H
#define OPENMC_CELL_H

#include <cstdint>
#include <functional> // for hash
#include <limits>
#include <string>
#include <unordered_map>
#include <unordered_set>

#include "hdf5.h"
#include "pugixml.hpp"

#include "openmc/bounding_box.h"
#include "openmc/constants.h"
#include "openmc/memory.h" // for unique_ptr
#include "openmc/neighbor_list.h"
#include "openmc/position.h"
#include "openmc/surface.h"
#include "openmc/universe.h"
#include "openmc/vector.h"

namespace openmc {

//==============================================================================
// Constants
//==============================================================================

enum class Fill { MATERIAL, UNIVERSE, LATTICE };

constexpr int32_t OP_LEFT_PAREN {std::numeric_limits<int32_t>::max()};
constexpr int32_t OP_RIGHT_PAREN {std::numeric_limits<int32_t>::max() - 1};
constexpr int32_t OP_COMPLEMENT {std::numeric_limits<int32_t>::max() - 2};
constexpr int32_t OP_INTERSECTION {std::numeric_limits<int32_t>::max() - 3};
constexpr int32_t OP_UNION {std::numeric_limits<int32_t>::max() - 4};

//==============================================================================
// Global variables
//==============================================================================

class Cell;
class GeometryState;
class ParentCell;
class CellInstance;
class Universe;
class UniversePartitioner;

namespace model {
extern std::unordered_map<int32_t, int32_t> cell_map;
extern vector<unique_ptr<Cell>> cells;

} // namespace model

//==============================================================================

class Region {
public:
  //----------------------------------------------------------------------------
  // Constructors
  Region() {}
  explicit Region(std::string region_spec, int32_t cell_id);

  //----------------------------------------------------------------------------
  // Methods

  //! \brief Determine if a cell contains the particle at a given location.
  //!
  //! The bounds of the cell are determined by a logical expression involving
  //! surface half-spaces, stored as an expression tree of intersections and
  //! unions with half-spaces as leaves.
  //!
  //! The function is split into two cases, one for simple cells (those
  //! involving only the intersection of half-spaces) and one for complex cells.
  //! Both cases use short circuiting.
  //! \param r The 3D Cartesian coordinate to check.
  //! \param u A direction used to "break ties" the coordinates are very
  //!   close to a surface.
  //! \param on_surface The signed index of a surface that the coordinate is
  //!   known to be on.  This index takes precedence over surface sense
  //!   calculations.
  bool contains(Position r, Direction u, int32_t on_surface) const;

  //! Find the oncoming boundary of this cell.
  //! \param p Particle whose scratch space is used for complex regions, or
  //!   nullptr outside of transport
  std::pair<double, int32_t> distance(Position r, Direction u,
    int32_t on_surface, GeometryState* p = nullptr) const;

  //! Get the BoundingBox for this cell.
  BoundingBox bounding_box() const;

  //! Get the CSG expression as a string
  std::string str() const;

  //! Get a vector containing all the half-spaces in the region expression
  vector<int32_t> surfaces() const;

  //! Get the number of half-spaces in the region expression
  int n_surfaces() const;

  //----------------------------------------------------------------------------
  // Accessors

  //! Get Boolean of if the cell is simple or not
  bool is_simple() const { return !complex_; }

private:
  //----------------------------------------------------------------------------
  // Types

  //! Node of the region expression tree. Nodes are stored in pre-order, so
  //! the children of an operator node follow it, and the subtree of a node
  //! ends just before index end. Children of an operator node are never
  //! operator nodes of the same type.
  struct Node {
    enum class Type : int8_t { HALFSPACE, INTERSECTION, UNION };
    Type type;
    int32_t halfspace; //!< Signed position + 1 in surfaces (HALFSPACE nodes)
    int32_t end;       //!< Index one past the last node of the subtree
    int32_t parent;    //!< Index of the parent node (-1 for the root)
  };

  //----------------------------------------------------------------------------
  // Private Methods

  //! Determine if a particle is inside the cell for a simple cell (only
  //! intersection operators)
  bool contains_simple(Position r, Direction u, int32_t on_surface) const;

  //! Determine if a particle is inside the cell for a complex cell.
  //!
  //! Evaluates the expression tree, skipping the remaining children of an
  //! operator node as soon as its value is known.
  bool contains_complex(Position r, Direction u, int32_t on_surface) const;

  //! Evaluate the expression tree of a complex region
  //!
  //! Operator nodes are evaluated with short circuiting, skipping their
  //! remaining children as soon as one of them determines their value.
  //! \param in_halfspace Callable returning whether the point is in the
  //!   half-space of a HALFSPACE node given its halfspace value
  template<typename F>
  bool evaluate(F&& in_halfspace) const;

  //! Signed surface index + 1 of a half-space of the expression tree
  int32_t surface_token(int32_t halfspace) const
  {
    int32_t i_surf = complex_->surfaces[std::abs(halfspace) - 1];
    return halfspace > 0 ? i_surf : -i_surf;
  }

  //! Find the nearest intersection with any surface in the region expression.
  std::pair<double, int32_t> distance_to_nearest_surface(Position r,
    Direction u, int32_t on_surface, bool ignore_coincident_surfaces) const;

  //! Find the oncoming boundary of this cell for a complex cell.
  std::pair<double, int32_t> distance_complex(
    Position r, Direction u, int32_t on_surface, GeometryState* p) const;

  //----------------------------------------------------------------------------
  // Private Data

  //! Signed surface indices + 1 of the half-spaces in the region expression,
  //! in order. A simple region is the intersection of these half-spaces.
  vector<int32_t> halfspaces_;

  //! Data needed only by complex regions, kept out of line so that regions,
  //! and the cells holding them, stay small for simple cells
  struct Complex {
    vector<Node> nodes; //!< Expression tree in pre-order
    //! Distinct surface indices + 1 of the half-spaces, in order of first
    //! appearance
    vector<int32_t> surfaces;
  };

  //! Data of a complex region (null for a simple region)
  unique_ptr<Complex> complex_;
};

//==============================================================================
// XML parsing helpers for <cell> nodes
//==============================================================================

//! Parse material IDs from a <cell> XML node.
//! \param node XML node containing a "material" attribute or child element
//! \param cell_id Cell ID used in error messages
//! \return Vector of material IDs (MATERIAL_VOID for "void")
vector<int32_t> parse_cell_material_xml(pugi::xml_node node, int32_t cell_id);

//! Parse temperatures in [K] from a <cell> XML node.
//! Validates that all values are non-negative and the list is non-empty.
//! \param node XML node containing a "temperature" attribute or child element
//! \param cell_id Cell ID used in error messages
//! \return Vector of temperatures in [K]
vector<double> parse_cell_temperature_xml(pugi::xml_node node, int32_t cell_id);

//! Parse densities in [g/cm³] from a <cell> XML node.
//! Validates that all values are positive and the list is non-empty.
//! \param node XML node containing a "density" attribute or child element
//! \param cell_id Cell ID used in error messages
//! \return Vector of densities in [g/cm³]
vector<double> parse_cell_density_xml(pugi::xml_node node, int32_t cell_id);

//==============================================================================

class Cell {
public:
  //----------------------------------------------------------------------------
  // Constructors, destructors, factory functions

  explicit Cell(pugi::xml_node cell_node);
  Cell() {};
  virtual ~Cell() = default;

  //----------------------------------------------------------------------------
  // Methods

  //! \brief Determine if a cell contains the particle at a given location.
  //!
  //! The bounds of the cell are detemined by a logical expression involving
  //! surface half-spaces. At initialization, the expression was converted
  //! to RPN notation.
  //!
  //! The function is split into two cases, one for simple cells (those
  //! involving only the intersection of half-spaces) and one for complex cells.
  //! Simple cells can be evaluated with short circuit evaluation, i.e., as soon
  //! as we know that one half-space is not satisfied, we can exit. This
  //! provides a performance benefit for the common case. In
  //! contains_complex, we evaluate the RPN expression using a stack, similar to
  //! how a RPN calculator would work.
  //! \param r The 3D Cartesian coordinate to check.
  //! \param u A direction used to "break ties" the coordinates are very
  //!   close to a surface.
  //! \param on_surface The signed index of a surface that the coordinate is
  //!   known to be on.  This index takes precedence over surface sense
  //!   calculations.
  virtual bool contains(Position r, Direction u, int32_t on_surface) const = 0;

  //! Find the oncoming boundary of this cell.
  virtual std::pair<double, int32_t> distance(
    Position r, Direction u, int32_t on_surface, GeometryState* p) const = 0;

  //! Write all information needed to reconstruct the cell to an HDF5 group.
  //! \param group_id An HDF5 group id.
  void to_hdf5(hid_t group_id) const;

  virtual void to_hdf5_inner(hid_t group_id) const = 0;

  //! Export physical properties to HDF5
  //! \param[in] group  HDF5 group to read from
  void export_properties_hdf5(hid_t group) const;

  //! Import physical properties from HDF5
  //! \param[in] group  HDF5 group to write to
  void import_properties_hdf5(hid_t group);

  //! Get the BoundingBox for this cell.
  virtual BoundingBox bounding_box() const = 0;

  //! Get a vector of surfaces in the cell
  virtual vector<int32_t> surfaces() const { return vector<int32_t>(); }

  //! Get the number of surfaces in the cell
  virtual int n_surfaces() const { return 0; }

  //! Check if the cell region expression is simple
  virtual bool is_simple() const { return true; }

  //----------------------------------------------------------------------------
  // Accessors

  //! Get the temperature of a cell instance
  //! \param[in] instance Instance index. If -1 is given, the temperature for
  //!   the first instance is returned.
  //! \return Temperature in [K]
  double temperature(int32_t instance = -1) const;

  //! Get the density multiplier of a cell instance
  //! \param[in] instance Instance index. If -1 is given, the density multiplier
  //! for the first instance is returned.
  //! \return Density multiplier
  double density_mult(int32_t instance = -1) const;

  //! Get the density of a cell instance in g/cm3
  //! \param[in] instance Instance index. If -1 is given, the density
  //! for the first instance is returned.
  //! \return Density in [g/cm3]
  double density(int32_t instance = -1) const;

  //! Set the temperature of a cell instance
  //! \param[in] T Temperature in [K]
  //! \param[in] instance Instance index. If -1 is given, the temperature for
  //!   all instances is set.
  //! \param[in] set_contained If this cell is not filled with a material,
  //!   collect all contained cells with material fills and set their
  //!   temperatures.
  void set_temperature(
    double T, int32_t instance = -1, bool set_contained = false);

  //! Set the density of a cell instance
  //! \param[in] density Density [g/cm3]
  //! \param[in] instance Instance index. If -1 is given, the density
  //!   for all instances is set.
  //! \param[in] set_contained If this cell is not filled with a material,
  //!   collect all contained cells with material fills and set their
  //!   densities.
  void set_density(
    double density, int32_t instance = -1, bool set_contained = false);

  int32_t n_instances() const;

  //! Set the rotation matrix of a cell instance
  //! \param[in] rot The rotation matrix of length 3 or 9
  void set_rotation(const vector<double>& rot);

  //! Get the name of a cell
  //! \return Cell name
  const std::string& name() const { return name_; };

  //! Set the temperature of a cell instance
  //! \param[in] name Cell name
  void set_name(const std::string& name) { name_ = name; };

  //! Get all cell instances contained by this cell
  //! \param[in] instance Instance of the cell for which to get contained cells
  //! (default instance is zero)
  //! \param[in] hint positional hint for determining the parent cells
  //! \return Map with cell indexes as keys and
  //! instances as values
  std::unordered_map<int32_t, vector<int32_t>> get_contained_cells(
    int32_t instance = 0, Position* hint = nullptr) const;

  //! Determine the material index corresponding to a specific cell instance,
  //! taking into account presence of distribcell material
  //! \param[in] instance of the cell
  //! \return material index
  int32_t material(int32_t instance) const
  {
    // If distributed materials are used, then each instance has its own
    // material definition. If distributed materials are not used, then
    // all instances used the same material stored at material_[0]. The
    // presence of distributed materials is inferred from the size of
    // the material_ vector being greater than one.
    if (material_.size() > 1) {
      return material_[instance];
    } else {
      return material_[0];
    }
  }

  //! Determine the temperature index corresponding to a specific cell instance,
  //! taking into account presence of distribcell temperature
  //! \param[in] instance of the cell
  //! \return temperature index
  double sqrtkT(int32_t instance) const
  {
    // If distributed materials are used, then each instance has its own
    // temperature definition. If distributed materials are not used, then
    // all instances used the same temperature stored at sqrtkT_[0]. The
    // presence of distributed materials is inferred from the size of
    // the sqrtkT_ vector being greater than one.
    if (sqrtkT_.size() > 1) {
      return sqrtkT_[instance];
    } else {
      return sqrtkT_[0];
    }
  }

protected:
  //! Determine the path to this cell instance in the geometry hierarchy
  //! \param[in] instance of the cell to find parent cells for
  //! \param[in] r position used to do a fast search for parent cells
  //! \return parent cells
  vector<ParentCell> find_parent_cells(
    int32_t instance, const Position& r) const;

  //! Determine the path to this cell instance in the geometry hierarchy
  //! \param[in] instance of the cell to find parent cells for
  //! \param[in] p particle used to do a fast search for parent cells
  //! \return parent cells
  vector<ParentCell> find_parent_cells(
    int32_t instance, GeometryState& p) const;

  //! Determine the path to this cell instance in the geometry hierarchy
  //! \param[in] instance of the cell to find parent cells for
  //! \return parent cells
  vector<ParentCell> exhaustive_find_parent_cells(int32_t instance) const;

  //! Inner function for retrieving contained cells
  void get_contained_cells_inner(
    std::unordered_map<int32_t, vector<int32_t>>& contained_cells,
    vector<ParentCell>& parent_cells) const;

public:
  //----------------------------------------------------------------------------
  // Data members

  int32_t id_;       //!< Unique ID
  std::string name_; //!< User-defined name
  Fill type_;        //!< Material, universe, or lattice
  int32_t universe_; //!< Universe # this cell is in
  int32_t fill_;     //!< Universe # filling this cell

  //! \brief Index corresponding to this cell in distribcell arrays
  int distribcell_index_ {C_NONE};

  //! \brief Material(s) within this cell.
  //!
  //! May be multiple materials for distribcell.
  vector<int32_t> material_;

  //! \brief Temperature(s) within this cell.
  //!
  //! The stored values are actually sqrt(k_Boltzmann * T) for each temperature
  //! T. The units are sqrt(eV).
  vector<double> sqrtkT_;

  //! \brief Unitless density multiplier(s) within this cell.
  vector<double> density_mult_;

  //! \brief Neighboring cells in the same universe.
  NeighborList neighbors_;

  Position translation_ {0, 0, 0}; //!< Translation vector for filled universe

  //! \brief Rotational tranfsormation of the filled universe.
  //
  //! The vector is empty if there is no rotation. Otherwise, the first 9 values
  //! give the rotation matrix in row-major order. When the user specifies
  //! rotation angles about the x-, y- and z- axes in degrees, these values are
  //! also present at the end of the vector, making it of length 12.
  vector<double> rotation_;

  vector<int32_t> offset_; //!< Distribcell offset table

  // Right now, either CSG or DAGMC cells are used.
  virtual GeometryType geom_type() const = 0;
};

struct CellInstanceItem {
  int32_t index {-1};    //! Index into global cells array
  int lattice_indx {-1}; //! Flat index value of the lattice cell
};

//==============================================================================

class CSGCell : public Cell {
public:
  //----------------------------------------------------------------------------
  // Constructors
  CSGCell() = default;
  explicit CSGCell(pugi::xml_node cell_node);

  //----------------------------------------------------------------------------
  // Methods
  vector<int32_t> surfaces() const override { return region_.surfaces(); }

  int n_surfaces() const override { return region_.n_surfaces(); }

  std::pair<double, int32_t> distance(Position r, Direction u,
    int32_t on_surface, GeometryState* p) const override
  {
    return region_.distance(r, u, on_surface, p);
  }

  bool contains(Position r, Direction u, int32_t on_surface) const override
  {
    return region_.contains(r, u, on_surface);
  }

  BoundingBox bounding_box() const override { return region_.bounding_box(); }

  void to_hdf5_inner(hid_t group_id) const override;

  bool is_simple() const override { return region_.is_simple(); }

  virtual GeometryType geom_type() const override { return GeometryType::CSG; }

private:
  Region region_;
};

//==============================================================================
//! Define an instance of a particular cell
//==============================================================================

//!  Stores information used to identify a unique cell in the model
struct CellInstance {
  //! Check for equality
  bool operator==(const CellInstance& other) const
  {
    return index_cell == other.index_cell && instance == other.instance;
  }

  int64_t index_cell;
  int64_t instance;
};

//! Structure necessary for inserting CellInstance into hashed STL data
//! structures
struct CellInstanceHash {
  std::size_t operator()(const CellInstance& k) const
  {
    return 4096 * k.index_cell + k.instance;
  }
};

//==============================================================================
// Non-member functions
//==============================================================================

void read_cells(pugi::xml_node node);

//!  Add cells to universes
void populate_universes();

} // namespace openmc
#endif // OPENMC_CELL_H
