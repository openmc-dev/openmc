#ifndef OPENMC_MATERIAL_H
#define OPENMC_MATERIAL_H

#include <memory> // for unique_ptr
#include <unordered_map>
#include <vector>

#include "pugixml.hpp"

#include "openmc/bremsstrahlung.h"
#include "openmc/particle.h"

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

class Material;

namespace model {

extern std::vector<Material*> materials;
extern std::unordered_map<int32_t, int32_t> material_map;

} // namespace model

//==============================================================================
//! A substance with constituent nuclides and thermal scattering data
//==============================================================================

class Material
{
public:
  int32_t id_; //!< Unique ID
  double volume_ {-1.0}; //!< Volume in [cm^3]
  bool fissionable {false}; //!< Does this material contain fissionable nuclides

  //! \brief Default temperature for cells containing this material.
  //!
  //! A negative value indicates no default temperature was specified.
  double temperature_ {-1};

  std::unique_ptr<Bremsstrahlung> ttb_;

  Material() {};

  explicit Material(pugi::xml_node material_node);

  //! Initialize bremsstrahlung data
  void init_bremsstrahlung();
};

//==============================================================================
// Fortran compatibility
//==============================================================================

extern "C" int* material_element(int i_material);
extern "C" bool material_isotropic(int i_material, int i_nuc_mat);

} // namespace openmc
#endif // OPENMC_MATERIAL_H
