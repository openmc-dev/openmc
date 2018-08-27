#include "openmc/material.h"

#include <string>
#include <sstream>

#include "openmc/error.h"
#include "openmc/xml_interface.h"


namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

std::vector<Material*> global_materials;
std::unordered_map<int32_t, int32_t> material_map;

//==============================================================================
// Material implementation
//==============================================================================

Material::Material(pugi::xml_node node)
{
  if (check_for_node(node, "id")) {
    id = std::stoi(get_node_value(node, "id"));
  } else {
    fatal_error("Must specify id of material in materials XML file.");
  }

  if (check_for_node(node, "volume")) {
    volume_ = std::stod(get_node_value(node, "volume"));
  }
}

//==============================================================================
// Non-method functions
//==============================================================================

extern "C" void
read_materials(pugi::xml_node* node)
{
  // Loop over XML material elements and populate the array.
  for (pugi::xml_node material_node : node->children("material")) {
    global_materials.push_back(new Material(material_node));
  }
  global_materials.shrink_to_fit();

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
// C API
//==============================================================================

extern "C" int
openmc_material_get_volume(int32_t index, double* volume)
{
  if (index >= 1 && index <= global_materials.size()) {
    Material* m = global_materials[index - 1];
    if (m->volume_ >= 0.0) {
      *volume = m->volume_;
      return 0;
    } else {
      std::stringstream msg;
      msg << "Volume for material with ID=" << m->id << " not set.";
      set_errmsg(msg);
      return OPENMC_E_UNASSIGNED;
    }
  } else {
    set_errmsg("Index in materials array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }
}

extern "C" int
openmc_material_set_volume(int32_t index, double volume)
{
  if (index >= 1 && index <= global_materials.size()) {
    Material* m = global_materials[index - 1];
    if (volume >= 0.0) {
      m->volume_ = volume;
      return 0;
    } else {
      set_errmsg("Volume must be non-negative");
      return OPENMC_E_INVALID_ARGUMENT;
    }
  } else {
    set_errmsg("Index in materials array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }
}

//==============================================================================
// Fortran compatibility functions
//==============================================================================

extern "C" {
  Material* material_pointer(int32_t indx) {return global_materials[indx];}

  int32_t material_id(Material* mat) {return mat->id;}

  void material_set_id(Material* mat, int32_t id, int32_t index)
  {
    mat->id = id;
    //TODO: off-by-one
    material_map[id] = index - 1;
  }

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
