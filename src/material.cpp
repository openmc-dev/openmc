#include "material.h"

#include "error.h"
#include "xml_interface.h"

//TODO: remove this include
#include <iostream>


namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

std::vector<Material*> global_materials;

//==============================================================================
// Material implementation
//==============================================================================

Material::Material(pugi::xml_node material_node)
{
  if (check_for_node(material_node, "id")) {
    id = stoi(get_node_value(material_node, "id"));
  } else {
    fatal_error("Must specify id of material in materials XML file.");
  }

  std::cout << "Reading material with id " << id << std::endl;
}

//==============================================================================
// Non-method functions
//==============================================================================

extern "C" void
read_materials(pugi::xml_node* node)
{
  // Loop over XML material elements and populate the array.
  for (pugi::xml_node material_node: node->children("material")) {
    global_materials.push_back(new Material(material_node));
  }
}

//==============================================================================
// Fortran compatibility functions
//==============================================================================

extern "C" {
  Material* material_pointer(int32_t indx) {return global_materials[indx];}

  int32_t material_id(Material* mat) {return mat->id;}

  void material_set_id(Material* mat, int32_t id) {mat->id = id;}

  void extend_materials_c(int32_t n)
  {
    global_materials.reserve(global_materials.size() + n);
    for (int32_t i = 0; i < n; i++) {
      global_materials.push_back(new Material());
    }
  }
}

} // namespace openmc
