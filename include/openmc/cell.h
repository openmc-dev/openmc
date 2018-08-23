#ifndef OPENMC_CELL_H
#define OPENMC_CELL_H

#include <cstdint>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

#include "hdf5.h"
#include "pugixml.hpp"

#include "openmc/position.h"


namespace openmc {

//==============================================================================
// Constants
//==============================================================================

// TODO: Convert to enum
extern "C" int FILL_MATERIAL;
extern "C" int FILL_UNIVERSE;
extern "C" int FILL_LATTICE;

//==============================================================================
// Global variables
//==============================================================================

extern "C" int32_t n_cells;

class Cell;
extern std::vector<Cell*> global_cells;
extern std::unordered_map<int32_t, int32_t> cell_map;

class Universe;
extern std::vector<Universe*> global_universes;
extern std::unordered_map<int32_t, int32_t> universe_map;

//==============================================================================
//! A geometry primitive that fills all space and contains cells.
//==============================================================================

class Universe
{
public:
  int32_t id;                  //!< Unique ID
  std::vector<int32_t> cells;  //!< Cells within this universe
  //double x0, y0, z0;           //!< Translation coordinates.
};

//==============================================================================
//! A geometry primitive that links surfaces, universes, and materials
//==============================================================================

class Cell
{
public:
  int32_t id;                //!< Unique ID
  std::string name;          //!< User-defined name
  int type;                  //!< Material, universe, or lattice
  int32_t universe;          //!< Universe # this cell is in
  int32_t fill;              //!< Universe # filling this cell
  int32_t n_instances{0};    //!< Number of instances of this cell

  //! \brief Material(s) within this cell.
  //!
  //! May be multiple materials for distribcell.  C_NONE signifies a universe.
  std::vector<int32_t> material;

  //! Definition of spatial region as Boolean expression of half-spaces
  std::vector<std::int32_t> region;
  //! Reverse Polish notation for region expression
  std::vector<std::int32_t> rpn;
  bool simple;  //!< Does the region contain only intersections?

  std::vector<int32_t> offset;  //!< Distribcell offset table

  Cell() {};

  explicit Cell(pugi::xml_node cell_node);

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
  distance(Position r, Direction u, int32_t on_surface) const;

  //! \brief Write cell information to an HDF5 group.
  //! \param group_id An HDF5 group id.
  void to_hdf5(hid_t group_id) const;

protected:
  bool contains_simple(Position r, Direction u, int32_t on_surface) const;
  bool contains_complex(Position r, Direction u, int32_t on_surface) const;
};

} // namespace openmc
#endif // OPENMC_CELL_H
