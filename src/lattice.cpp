#include "lattice.h"

#include <sstream>
#include <vector>

#include "error.h"
#include "hdf5_interface.h"
#include "xml_interface.h"

//TODO: remove this include
#include <iostream>


namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

// Braces force n_lattices to be defined here, not just declared.
//extern "C" {int32_t n_lattices {0};}

//Lattice **lattices_c;
std::vector<Lattice*> lattices_c;

std::map<int32_t, int32_t> lattice_dict;

//==============================================================================
// Lattice implementation
//==============================================================================

Lattice::Lattice(pugi::xml_node lat_node)
{
  if (check_for_node(lat_node, "id")) {
    id = stoi(get_node_value(lat_node, "id"));
  } else {
    fatal_error("Must specify id of lattice in geometry XML file.");
  }

  if (check_for_node(lat_node, "name")) {
    name = get_node_value(lat_node, "name");
  }
}

void
Lattice::to_hdf5(hid_t lat_group) const
{
  if (!name.empty()) {
    write_string(lat_group, "name", name);
  }
}

//==============================================================================
// RectLattice implementation
//==============================================================================

RectLattice::RectLattice(pugi::xml_node lat_node)
  : Lattice {lat_node}
{
}

//==============================================================================
// HexLattice implementation
//==============================================================================

HexLattice::HexLattice(pugi::xml_node lat_node)
  : Lattice {lat_node}
{
}

//==============================================================================

extern "C" void
read_lattices(pugi::xml_node *node)
{
  for (pugi::xml_node lat_node: node->children("lattice")) {
    lattices_c.push_back(new RectLattice(lat_node));
  }
  for (pugi::xml_node lat_node: node->children("hex_lattice")) {
    lattices_c.push_back(new HexLattice(lat_node));
  }

  // Fill the lattice dictionary.
  for (int i_lat = 0; i_lat < lattices_c.size(); i_lat++) {
    int id = lattices_c[i_lat]->id;
    auto in_dict = lattice_dict.find(id);
    if (in_dict == lattice_dict.end()) {
      lattice_dict[id] = i_lat;
    } else {
      std::stringstream err_msg;
      err_msg << "Two or more lattices use the same unique ID: " << id;
      fatal_error(err_msg);
    }
  }
}

//==============================================================================
// Fortran compatibility functions
//==============================================================================

extern "C" {
  Lattice* lattice_pointer(int lat_ind) {return lattices_c[lat_ind];}

  int32_t lattice_id(Lattice *lat) {return lat->id;}

  void lattice_to_hdf5(Lattice *lat, hid_t group) {lat->to_hdf5(group);}
}

} // namespace openmc
