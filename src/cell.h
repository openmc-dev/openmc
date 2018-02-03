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
// Module constants
//==============================================================================

//==============================================================================
// Global variables
//==============================================================================

// Braces force n_cells to be defined here, not just declared.
extern "C" {int32_t n_cells {0};}

class Cell;
//Cell *cells_c;
std::vector<Cell> cells_c;

std::map<int, int> cell_dict;

//==============================================================================
//! A geometry primitive that links surfaces, universes, and materials
//==============================================================================

class Cell 
{
public:
  int32_t id;                //!< Unique ID
  std::string name{""};      //!< User-defined name

  explicit Cell(pugi::xml_node cell_node);

  //! Write all information needed to reconstruct the cell to an HDF5 group.
  //! @param group_id An HDF5 group id.
  void to_hdf5(hid_t group_id) const;
};

//==============================================================================
// Fortran compatibility functions
//==============================================================================

extern "C" void read_cells(pugi::xml_node *node);

extern "C" Cell* cell_pointer(int cell_ind) {return &cells_c[cell_ind];}

extern "C" int32_t cell_id(Cell *c) {return c->id;}

extern "C" void cell_set_id(Cell *c, int32_t id) {c->id = id;}

//extern "C" void free_memory_cells_c()
//{
//  delete cells_c;
//  cells_c = nullptr;
//  n_cells = 0;
//  cell_dict.clear();
//}

} // namespace openmc
#endif // CELL_H
