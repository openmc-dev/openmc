#ifndef OPENMC_MATERIAL_H
#define OPENMC_MATERIAL_H

#include <memory> // for unique_ptr
#include <string>
#include <unordered_map>
#include <vector>

#include "pugixml.hpp"
#include "xtensor/xtensor.hpp"

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
  // Constructors
  Material() {};
  explicit Material(pugi::xml_node material_node);

  // Methods
  void calculate_xs(const Particle& p) const;

  //! Initialize bremsstrahlung data
  void init_bremsstrahlung();

  // Data
  int32_t id_; //!< Unique ID
  std::string name_; //!< Name of material
  xt::xtensor<int, 1> nuclide_; //!< Indices in nuclides vector
  xt::xtensor<int, 1> element_; //!< Indices in elements vector
  xt::xtensor<double, 1> atom_density_; //!< Nuclide atom density in [atom/b-cm]
  double density_; //!< Total atom density in [atom/b-cm]
  double volume_ {-1.0}; //!< Volume in [cm^3]
  bool fissionable_ {false}; //!< Does this material contain fissionable nuclides
  bool depletable_ {false}; //!< Is the material depletable?

  //! \brief Default temperature for cells containing this material.
  //!
  //! A negative value indicates no default temperature was specified.
  double temperature_ {-1};

  std::unique_ptr<Bremsstrahlung> ttb_;

private:
  void calculate_neutron_xs(const Particle& p) const;
  void calculate_photon_xs(const Particle& p) const;
};

//==============================================================================
// Fortran compatibility
//==============================================================================

extern "C" int* material_element(int i_material);
extern "C" bool material_isotropic(int i_material, int i_nuc_mat);

} // namespace openmc
#endif // OPENMC_MATERIAL_H
