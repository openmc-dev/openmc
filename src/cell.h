#ifndef CELL_H
#define CELL_H

#include <map>
#include <cstdint>
#include <string>
#include <vector>

#include "hdf5.h"
#include "pugixml/pugixml.hpp"


namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

// Braces force n_cells to be defined here, not just declared.
extern "C" {int32_t n_cells {0};}

class Cell;
std::vector<Cell*> cells_c;
std::map<int32_t, int32_t> cell_dict;

class Universe;
std::vector<Universe*> universes_c;
std::map<int32_t, int32_t> universe_dict;

//==============================================================================
//! A geometry primitive that fills all space and contains cells.
//==============================================================================

class Universe
{
public:
  int32_t id;                  //! Unique ID
  int32_t type;
  std::vector<Cell*> cells;    //! Cells within this universe
  double x0, y0, z0;           //! Translation coordinates.
};

//==============================================================================
//! A geometry primitive that links surfaces, universes, and materials
//==============================================================================

class Cell 
{
public:
  int32_t id;                //!< Unique ID
  std::string name{""};      //!< User-defined name
  int32_t universe;          //!< Universe # this cell is in
  int32_t fill;              //!< Universe # filling this cell

  //! Definition of spatial region as Boolean expression of half-spaces
  std::vector<std::int32_t> region;
  //! Reverse Polish notation for region expression
  std::vector<std::int32_t> rpn;
  bool simple;  //!< Does the region contain only intersections?

  Cell() {};

  explicit Cell(pugi::xml_node cell_node);

  //! Determine if a cell contains the particle at a given location.
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
  //! contains_complex, we evaluate the RPN expression using a stack, similar to  //! how a RPN calculator would work.
  //! @param xyz[3] The 3D Cartesian coordinate to check.
  //! @param uvw[3] A direction used to "break ties" the coordinates are very
  //!   close to a surface.
  //! @param on_surface The signed index of a surface that the coordinate is
  //!   known to be on.  This index takes precedence over surface sense
  //!   calculations.
  bool
  contains(const double xyz[3], const double uvw[3], int32_t on_surface) const;

  std::pair<double, int32_t>
  distance(const double xyz[3], const double uvw[3], int32_t on_surface) const;

  //! Write all information needed to reconstruct the cell to an HDF5 group.
  //! @param group_id An HDF5 group id.
  void to_hdf5(hid_t group_id) const;

protected:
  bool contains_simple(const double xyz[3], const double uvw[3],
                       int32_t on_surface) const;
  bool contains_complex(const double xyz[3], const double uvw[3],
                        int32_t on_surface) const;
};

} // namespace openmc
#endif // CELL_H
