#include "cell.h"

#include <sstream>

#include "error.h"
#include "hdf5_interface.h"
#include "xml_interface.h"

//TODO: remove this include
#include <iostream>


namespace openmc {

//==============================================================================
// Cell implementation
//==============================================================================

Cell::Cell(pugi::xml_node cell_node)
{
  if (check_for_node(cell_node, "id")) {
    id = stoi(get_node_value(cell_node, "id"));
  } else {
    fatal_error("Must specify id of cell in geometry XML file.");
  }

  //TODO: don't automatically lowercase cell and surface names
  if (check_for_node(cell_node, "name")) {
    name = get_node_value(cell_node, "name");
  }
}

void
Cell::to_hdf5(hid_t group_id) const
{
/*
  std::string group_name {"surface "};
  group_name += std::to_string(id);

  hid_t surf_group = create_group(group_id, group_name);

  switch(bc) {
    case BC_TRANSMIT :
      write_string(surf_group, "boundary_type", "transmission");
      break;
    case BC_VACUUM :
      write_string(surf_group, "boundary_type", "vacuum");
      break;
    case BC_REFLECT :
      write_string(surf_group, "boundary_type", "reflective");
      break;
    case BC_PERIODIC :
      write_string(surf_group, "boundary_type", "periodic");
      break;
  }

  if (!name.empty()) {
    write_string(surf_group, "name", name);
  }

  to_hdf5_inner(surf_group);

  close_group(surf_group);
*/
}

//==============================================================================

extern "C" void
read_cells(pugi::xml_node *node)
{
  // Count the number of cells.
  for (pugi::xml_node cell_node: node->children("cell")) {n_cells++;}
  if (n_cells == 0) {
    fatal_error("No cells found in geometry.xml!");
  }

  // Allocate the vector of Cells.
  cells_c.reserve(n_cells);

  // Loop over XML cell elements and populate the array.
  for (pugi::xml_node cell_node: node->children("cell")) {
    cells_c.push_back(Cell(cell_node));
  }

}

} // namespace openmc
