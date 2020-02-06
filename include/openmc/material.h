#ifndef OPENMC_MATERIAL_H
#define OPENMC_MATERIAL_H

#include <memory> // for unique_ptr
#include <string>
#include <unordered_map>
#include <vector>

#include <gsl/gsl>
#include <hdf5.h>
#include "pugixml.hpp"
#include "xtensor/xtensor.hpp"

#include "openmc/constants.h"
#include "openmc/bremsstrahlung.h"
#include "openmc/particle.h"

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

class Material;

namespace model {

extern std::unordered_map<int32_t, int32_t> material_map;
extern std::vector<std::unique_ptr<Material>> materials;

} // namespace model

//==============================================================================
//! A substance with constituent nuclides and thermal scattering data
//==============================================================================

class Material
{
public:
  //----------------------------------------------------------------------------
  // Types
  struct ThermalTable {
    int index_table; //!< Index of table in data::thermal_scatt
    int index_nuclide; //!< Index in nuclide_
    double fraction; //!< How often to use table
  };

  //----------------------------------------------------------------------------
  // Constructors, destructors, factory functions
  Material() {};
  explicit Material(pugi::xml_node material_node);
  ~Material();

  //----------------------------------------------------------------------------
  // Methods

  void calculate_xs(Particle& p) const;

  //! Assign thermal scattering tables to specific nuclides within the material
  //! so the code knows when to apply bound thermal scattering data
  void init_thermal();

  //! Set up mapping between global nuclides vector and indices in nuclide_
  void init_nuclide_index();

  //! Finalize the material, assigning tables, normalize density, etc.
  void finalize();

  //! Write material data to HDF5
  void to_hdf5(hid_t group) const;

  //! Add nuclide to the material
  //
  //! \param[in] nuclide Name of the nuclide
  //! \param[in] density Density of the nuclide in [atom/b-cm]
  void add_nuclide(const std::string& nuclide, double density);

  //! Set atom densities for the material
  //
  //! \param[in] name Name of each nuclide
  //! \param[in] density Density of each nuclide in [atom/b-cm]
  void set_densities(const std::vector<std::string>& name,
    const std::vector<double>& density);

  //----------------------------------------------------------------------------
  // Accessors

  //! Get density in [atom/b-cm]
  //! \return Density in [atom/b-cm]
  double density() const { return density_; }

  //! Get density in [g/cm^3]
  //! \return Density in [g/cm^3]
  double density_gpcc() const { return density_gpcc_; }

  //! Get name
  //! \return Material name
  const std::string& name() const { return name_; }

  //! Set name
  void set_name(const std::string& name) { name_ = name; }

  //! Set total density of the material
  //
  //! \param[in] density Density value
  //! \param[in] units Units of density
  void set_density(double density, gsl::cstring_span units);

  //! Get nuclides in material
  //! \return Indices into the global nuclides vector
  gsl::span<const int> nuclides() const { return {nuclide_.data(), nuclide_.size()}; }

  //! Get densities of each nuclide in material
  //! \return Densities in [atom/b-cm]
  gsl::span<const double> densities() const { return {atom_density_.data(), atom_density_.size()}; }

  //! Get ID of material
  //! \return ID of material
  int32_t id() const { return id_; }

  //! Assign a unique ID to the material
  //! \param[in] Unique ID to assign. A value of -1 indicates that an ID
  //!   should be automatically assigned.
  void set_id(int32_t id);

  //! Get whether material is fissionable
  //! \return Whether material is fissionable
  bool fissionable() const { return fissionable_; }

  //! Get volume of material
  //! \return Volume in [cm^3]
  double volume() const;

  //----------------------------------------------------------------------------
  // Data
  int32_t id_ {C_NONE}; //!< Unique ID
  std::string name_; //!< Name of material
  std::vector<int> nuclide_; //!< Indices in nuclides vector
  std::vector<int> element_; //!< Indices in elements vector
  xt::xtensor<double, 1> atom_density_; //!< Nuclide atom density in [atom/b-cm]
  double density_; //!< Total atom density in [atom/b-cm]
  double density_gpcc_; //!< Total atom density in [g/cm^3]
  double volume_ {-1.0}; //!< Volume in [cm^3]
  bool fissionable_ {false}; //!< Does this material contain fissionable nuclides
  bool depletable_ {false}; //!< Is the material depletable?
  std::vector<bool> p0_; //!< Indicate which nuclides are to be treated with iso-in-lab scattering

  // To improve performance of tallying, we store an array (direct address
  // table) that indicates for each nuclide in data::nuclides the index of the
  // corresponding nuclide in the nuclide_ vector. If it is not present in the
  // material, the entry is set to -1.
  std::vector<int> mat_nuclide_index_;

  // Thermal scattering tables
  std::vector<ThermalTable> thermal_tables_;

  //! \brief Default temperature for cells containing this material.
  //!
  //! A negative value indicates no default temperature was specified.
  double temperature_ {-1};

  std::unique_ptr<Bremsstrahlung> ttb_;

private:
  //----------------------------------------------------------------------------
  // Private methods

  //! Calculate the collision stopping power
  void collision_stopping_power(double* s_col, bool positron);

  //! Initialize bremsstrahlung data
  void init_bremsstrahlung();

  //! Normalize density
  void normalize_density();

  void calculate_neutron_xs(Particle& p) const;
  void calculate_photon_xs(Particle& p) const;

  //----------------------------------------------------------------------------
  // Private data members
  gsl::index index_;
};

//==============================================================================
// Non-member functions
//==============================================================================

//! Calculate Sternheimer adjustment factor
double sternheimer_adjustment(const std::vector<double>& f, const
  std::vector<double>& e_b_sq, double e_p_sq, double n_conduction, double
  log_I, double tol, int max_iter);

//! Calculate density effect correction
double density_effect(const std::vector<double>& f, const std::vector<double>&
  e_b_sq, double e_p_sq, double n_conduction, double rho, double E, double tol,
  int max_iter);

//! Read material data from materials.xml
void read_materials_xml();

void free_memory_material();

} // namespace openmc
#endif // OPENMC_MATERIAL_H
