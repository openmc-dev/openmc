#include "cell.h"

#include <sstream>
#include <string>

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

  if (check_for_node(cell_node, "universe")) {
    universe = stoi(get_node_value(cell_node, "universe"));
  } else {
    universe = 0;
  }
}

void
//Cell::to_hdf5(hid_t group_id) const
Cell::to_hdf5(hid_t cell_group) const
{
//  std::string group_name {"surface "};
//  group_name += std::to_string(id);
//
//  hid_t surf_group = create_group(group_id, group_name);

  if (!name.empty()) {
    write_string(cell_group, "name", name);
  }

  //TODO: Lookup universe id in universe_dict
  //write_int(cell_group, "universe", universe);

//  close_group(cell_group);
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
