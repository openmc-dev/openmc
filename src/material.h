#ifndef OPENMC_CELL_H
#define OPENMC_CELL_H

#include <vector>

#include "pugixml.hpp"


namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

class Material;
extern std::vector<Material*> global_materials;

//==============================================================================
//! A substance with constituent nuclides and thermal scattering data
//==============================================================================

class Material
{
public:
  int32_t id; //!< Unique ID

  Material() {};

  explicit Material(pugi::xml_node material_node);
};

} // namespace openmc
#endif // OPENMC_CELL_H
