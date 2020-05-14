#ifndef OPENMC_CELL_H
#define OPENMC_CELL_H

#include <cstdint>
#include <functional> // for hash
#include <limits>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include <gsl/gsl>
#include "hdf5.h"
#include "pugixml.hpp"
#include "dagmc.h"

#include "openmc/constants.h"
#include "openmc/neighbor_list.h"
#include "openmc/position.h"
#include "openmc/surface.h"

// Regression tests
/*
#define MATERIAL_SIZE 9 // CMFD NG
#define SQRTKT_SIZE 9 // CMFD NG
#define REGION_SIZE 32 // complex cell
#define RPN_SIZE 24 // complex cell
#define ROTATION_SIZE 12
#define OFFSET_SIZE 35
*/

// SMR
#define MATERIAL_SIZE 31680 // CMFD NG
#define SQRTKT_SIZE 31680 // CMFD NG
#define REGION_SIZE 32 // complex cell
#define RPN_SIZE 24 // complex cell
#define ROTATION_SIZE 12
#define OFFSET_SIZE 280

namespace openmc {

//==============================================================================
// Constants
//==============================================================================

enum class Fill {
  MATERIAL,
  UNIVERSE,
  LATTICE
};

// TODO: Convert to enum
constexpr int32_t OP_LEFT_PAREN   {std::numeric_limits<int32_t>::max()};
constexpr int32_t OP_RIGHT_PAREN  {std::numeric_limits<int32_t>::max() - 1};
constexpr int32_t OP_COMPLEMENT   {std::numeric_limits<int32_t>::max() - 2};
constexpr int32_t OP_INTERSECTION {std::numeric_limits<int32_t>::max() - 3};
constexpr int32_t OP_UNION        {std::numeric_limits<int32_t>::max() - 4};

//==============================================================================
// Global variables
//==============================================================================

class Cell;
class Universe;
class UniversePartitioner;

namespace model {
  //extern std::vector<std::unique_ptr<Cell>> cells;
  extern std::vector<Cell> cells;
  #pragma omp declare target
  extern Cell* device_cells;
  #pragma omp end declare target
  extern std::unordered_map<int32_t, int32_t> cell_map;

  //extern std::vector<std::unique_ptr<Universe>> universes;
  extern std::vector<Universe> universes;
  extern Universe* device_universes;
  extern std::unordered_map<int32_t, int32_t> universe_map;
} // namespace model

//==============================================================================
//! A geometry primitive that fills all space and contains cells.
//==============================================================================

class Universe
{
public:
  int32_t id_;                  //!< Unique ID
  std::vector<int32_t> cells_;  //!< Cells within this universe
  int32_t* device_cells_;  //!< Cells within this universe

  //! \brief Write universe information to an HDF5 group.
  //! \param group_id An HDF5 group id.
  void to_hdf5(hid_t group_id) const;
  
  void allocate_and_copy_to_device();

  BoundingBox bounding_box() const;

  std::unique_ptr<UniversePartitioner> partitioner_;
  UniversePartitioner* device_partitioner_{NULL};
};

//==============================================================================
//! A geometry primitive that links surfaces, universes, and materials
//==============================================================================

class Cell {
public:
  //----------------------------------------------------------------------------
  // Constructors, destructors, factory functions

  explicit Cell(pugi::xml_node cell_node);
  Cell() {};
  ~Cell() = default;

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
  bool
  contains(Position r, Direction u, int32_t on_surface) const;

  //! Find the oncoming boundary of this cell.
  std::pair<double, int32_t>
  distance(Position r, Direction u, int32_t on_surface, Particle* p) const;

  //! Write all information needed to reconstruct the cell to an HDF5 group.
  //! \param group_id An HDF5 group id.
  void to_hdf5(hid_t group_id) const;

  //! Get the BoundingBox for this cell.
  BoundingBox bounding_box() const;

  void allocate_on_device();
  void copy_to_device();

  //----------------------------------------------------------------------------
  // Accessors

  //! Get the temperature of a cell instance
  //! \param[in] instance Instance index. If -1 is given, the temperature for
  //!   the first instance is returned.
  //! \return Temperature in [K]
  double temperature(int32_t instance = -1) const;

  //! Set the temperature of a cell instance
  //! \param[in] T Temperature in [K]
  //! \param[in] instance Instance index. If -1 is given, the temperature for
  //!   all instances is set.
  void set_temperature(double T, int32_t instance = -1);

  //! Get the name of a cell
  //! \return Cell name
  const std::string& name() const { return name_; };

  //! Set the temperature of a cell instance
  //! \param[in] name Cell name
  void set_name(const std::string& name) { name_ = name; };

  //----------------------------------------------------------------------------
  // Data members

  int32_t id_;                //!< Unique ID
  std::string name_;          //!< User-defined name
  Fill type_;                  //!< Material, universe, or lattice
  int32_t universe_;          //!< Universe # this cell is in
  int32_t fill_;              //!< Universe # filling this cell
  int32_t n_instances_{0};    //!< Number of instances of this cell

  //! \brief Index corresponding to this cell in distribcell arrays
  int distribcell_index_{C_NONE};

  //! \brief Material(s) within this cell.
  //!
  //! May be multiple materials for distribcell.
  std::vector<int32_t> material_;
  int32_t* device_material_{NULL};

  //! \brief Temperature(s) within this cell.
  //!
  //! The stored values are actually sqrt(k_Boltzmann * T) for each temperature
  //! T. The units are sqrt(eV).
  std::vector<double> sqrtkT_;
  double* device_sqrtkT_{NULL};

  //! Definition of spatial region as Boolean expression of half-spaces
  std::vector<std::int32_t> region_;
  int32_t* device_region_{NULL};
  //! Reverse Polish notation for region expression
  std::vector<std::int32_t> rpn_;
  int32_t* device_rpn_{NULL};
  bool simple_;  //!< Does the region contain only intersections?

  //! \brief Neighboring cells in the same universe.
  NeighborList neighbors_;

  Position translation_ {0, 0, 0}; //!< Translation vector for filled universe

  //! \brief Rotational tranfsormation of the filled universe.
  //
  //! The vector is empty if there is no rotation. Otherwise, the first 9 values
  //! give the rotation matrix in row-major order. When the user specifies
  //! rotation angles about the x-, y- and z- axes in degrees, these values are
  //! also present at the end of the vector, making it of length 12.
  std::vector<double> rotation_;
  double* device_rotation_{NULL};

  std::vector<int32_t> offset_;  //!< Distribcell offset table
  int32_t* device_offset_{NULL};

protected:
  bool contains_simple(Position r, Direction u, int32_t on_surface) const;
  bool contains_complex(Position r, Direction u, int32_t on_surface) const;
  BoundingBox bounding_box_simple() const;
  static BoundingBox bounding_box_complex(std::vector<int32_t> rpn);

  //! Applies DeMorgan's laws to a section of the RPN
  //! \param start Starting point for token modification
  //! \param stop Stopping point for token modification
  static void apply_demorgan(std::vector<int32_t>::iterator start,
                             std::vector<int32_t>::iterator stop);

  //! Removes complement operators from the RPN
  //! \param rpn The rpn to remove complement operators from.
  static void remove_complement_ops(std::vector<int32_t>& rpn);

  //! Returns the beginning position of a parenthesis block (immediately before
  //! two surface tokens) in the RPN given a starting position at the end of
  //! that block (immediately after two surface tokens)
  //! \param start Starting position of the search
  //! \param rpn The rpn being searched
  static std::vector<int32_t>::iterator
  find_left_parenthesis(std::vector<int32_t>::iterator start,
                        const std::vector<int32_t>& rpn);
};

//==============================================================================

/*
class CSGCell : public Cell
{
public:
  CSGCell();

  explicit CSGCell(pugi::xml_node cell_node);

  bool
  contains(Position r, Direction u, int32_t on_surface) const;

  std::pair<double, int32_t>
  distance(Position r, Direction u, int32_t on_surface, Particle* p) const;

  void to_hdf5(hid_t group_id) const;

  BoundingBox bounding_box() const;

protected:
  bool contains_simple(Position r, Direction u, int32_t on_surface) const;
  bool contains_complex(Position r, Direction u, int32_t on_surface) const;
  BoundingBox bounding_box_simple() const;
  static BoundingBox bounding_box_complex(std::vector<int32_t> rpn);

  //! Applies DeMorgan's laws to a section of the RPN
  //! \param start Starting point for token modification
  //! \param stop Stopping point for token modification
  static void apply_demorgan(std::vector<int32_t>::iterator start,
                             std::vector<int32_t>::iterator stop);

  //! Removes complement operators from the RPN
  //! \param rpn The rpn to remove complement operators from.
  static void remove_complement_ops(std::vector<int32_t>& rpn);

  //! Returns the beginning position of a parenthesis block (immediately before
  //! two surface tokens) in the RPN given a starting position at the end of
  //! that block (immediately after two surface tokens)
  //! \param start Starting position of the search
  //! \param rpn The rpn being searched
  static std::vector<int32_t>::iterator
  find_left_parenthesis(std::vector<int32_t>::iterator start,
                        const std::vector<int32_t>& rpn);

};
                        */

//==============================================================================

#ifdef DAGMC
class DAGCell : public Cell
{
public:
  DAGCell();

  bool contains(Position r, Direction u, int32_t on_surface) const;

  std::pair<double, int32_t>
  distance(Position r, Direction u, int32_t on_surface, Particle* p) const;

  BoundingBox bounding_box() const;

  void to_hdf5(hid_t group_id) const;

  moab::DagMC* dagmc_ptr_; //!< Pointer to DagMC instance
  int32_t dag_index_;      //!< DagMC index of cell
};
#endif

//==============================================================================
//! Speeds up geometry searches by grouping cells in a search tree.
//
//! Currently this object only works with universes that are divided up by a
//! bunch of z-planes.  It could be generalized to other planes, cylinders,
//! and spheres.
//==============================================================================

class UniversePartitioner
{
public:
  explicit UniversePartitioner(const Universe& univ);

  //! Return the list of cells that could contain the given coordinates.
  //const std::vector<int32_t>& get_cells(Position r, Direction u) const;
  int32_t* get_cells(Position r, Direction u, int& ncells) const;

//private:
  //! A sorted vector of indices to surfaces that partition the universe
  std::vector<int32_t> surfs_;
  int32_t* device_surfs_{NULL};

  //! Vectors listing the indices of the cells that lie within each partition
  //
  //! There are n+1 partitions with n surfaces.  `partitions_.front()` gives the
  //! cells that lie on the negative side of `surfs_.front()`.
  //! `partitions_.back()` gives the cells that lie on the positive side of
  //! `surfs_.back()`.  Otherwise, `partitions_[i]` gives cells sandwiched
  //! between `surfs_[i-1]` and `surfs_[i]`.
  std::vector<std::vector<int32_t>> partitions_;
  int32_t** device_partitions_{NULL};
  int32_t* device_partitions_lengths_{NULL};
};

//==============================================================================
//! Define an instance of a particular cell
//==============================================================================

struct CellInstance {
  //! Check for equality
  bool operator==(const CellInstance& other) const
  { return index_cell == other.index_cell && instance == other.instance; }

  gsl::index index_cell;
  gsl::index instance;
};

struct CellInstanceHash {
  std::size_t operator()(const CellInstance& k) const
  {
    return 4096*k.index_cell + k.instance;
  }
};

//==============================================================================
// Non-member functions
//==============================================================================

void read_cells(pugi::xml_node node);

#ifdef DAGMC
int32_t next_cell(DAGCell* cur_cell, DAGSurface* surf_xed);
#endif

} // namespace openmc
#endif // OPENMC_CELL_H
