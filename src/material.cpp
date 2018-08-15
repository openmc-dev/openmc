#include "material.h"

#include <string>
#include <sstream>

#include "error.h"
#include "xml_interface.h"

//TODO: remove this include
#include <iostream>


namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

std::vector<Material*> global_materials;
std::unordered_map<int32_t, int32_t> material_map;

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

  // Populate the material map.
  for (int i = 0; i < global_materials.size(); i++) {
    int32_t mid = global_materials[i]->id;
    auto search = material_map.find(mid);
    if (search == material_map.end()) {
      material_map[mid] = i;
    } else {
      std::stringstream err_msg;
      err_msg << "Two or more materials use the same unique ID: " << mid;
      fatal_error(err_msg);
    }
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

  void free_memory_material_c()
  {
    for (Material *mat : global_materials) {delete mat;}
    global_materials.clear();
    material_map.clear();
  }
}

} // namespace openmc
