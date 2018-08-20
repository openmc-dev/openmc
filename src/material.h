#ifndef OPENMC_MATERIAL_H
#define OPENMC_MATERIAL_H

#include <unordered_map>
#include <vector>

#include "pugixml.hpp"


namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

class Material;
extern std::vector<Material*> global_materials;
extern std::unordered_map<int32_t, int32_t> material_map;

//==============================================================================
//! A substance with constituent nuclides and thermal scattering data
//==============================================================================

class Material
{
public:
  int32_t id; //!< Unique ID
  double volume_ {-1.0}; //!< Volume in [cm^3]

  Material() {};

  explicit Material(pugi::xml_node material_node);
};

} // namespace openmc
#endif // OPENMC_MATERIAL_H
