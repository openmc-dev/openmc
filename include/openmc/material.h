#ifndef OPENMC_MATERIAL_H
#define OPENMC_MATERIAL_H

#include <string>
#include <unordered_map>

#include "openmc/span.h"
#include "pugixml.hpp"
#include "xtensor/xtensor.hpp"
#include <hdf5.h>

#include "openmc/bremsstrahlung.h"
#include "openmc/constants.h"
#include "openmc/memory.h" // for unique_ptr
#include "openmc/ncrystal_interface.h"
#include "openmc/particle.h"
#include "openmc/vector.h"

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

class Material;
class Stochastic_Media;

namespace model {

extern std::unordered_map<int32_t, int32_t> material_map;
extern vector<unique_ptr<Material>> materials;

extern std::unordered_map<int32_t, int32_t> stochastic_media_map;
extern vector<unique_ptr<Stochastic_Media>> stochastic_media;

} // namespace model

//==============================================================================
//! A substance with constituent nuclides and thermal scattering data
//==============================================================================

class Material {
public:
  //----------------------------------------------------------------------------
  // Types
  struct ThermalTable {
    int index_table;   //!< Index of table in data::thermal_scatt
    int index_nuclide; //!< Index in nuclide_
    double fraction;   //!< How often to use table
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

  //! Export physical properties to HDF5
  //! \param[in] group  HDF5 group to write to
  void export_properties_hdf5(hid_t group) const;

  //! Import physical properties from HDF5
  //! \param[in] group  HDF5 group to read from
  void import_properties_hdf5(hid_t group);

  //! Add nuclide to the material
  //
  //! \param[in] nuclide Name of the nuclide
  //! \param[in] density Density of the nuclide in [atom/b-cm]
  void add_nuclide(const std::string& nuclide, double density);

  //! Set atom densities for the material
  //
  //! \param[in] name Name of each nuclide
  //! \param[in] density Density of each nuclide in [atom/b-cm]
  void set_densities(
    const vector<std::string>& name, const vector<double>& density);

  //! Clone the material by deep-copying all members, except for the ID,
  //  which will get auto-assigned to the next available ID. After creating
  //  the new material, it is added to openmc::model::materials.
  //! \return reference to the cloned material
  Material& clone();

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
  void set_density(double density, const std::string& units);

  //! Set temperature of the material
  void set_temperature(double temperature) { temperature_ = temperature; };

  //! Get nuclides in material
  //! \return Indices into the global nuclides vector
  span<const int> nuclides() const
  {
    return {nuclide_.data(), nuclide_.size()};
  }

  //! Get densities of each nuclide in material
  //! \return Densities in [atom/b-cm]
  span<const double> densities() const
  {
    return {atom_density_.data(), atom_density_.size()};
  }

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
  bool& fissionable() { return fissionable_; }

  //! Get volume of material
  //! \return Volume in [cm^3]
  double volume() const;

  //! Get temperature of material
  //! \return Temperature in [K]
  double temperature() const;

  //! Whether or not the material is depletable
  bool depletable() const { return depletable_; }
  bool& depletable() { return depletable_; }

  //! Get pointer to NCrystal material object
  //! \return Pointer to NCrystal material object
  const NCrystalMat& ncrystal_mat() const { return ncrystal_mat_; };

  //----------------------------------------------------------------------------
  // Data
  int32_t id_ {C_NONE};                 //!< Unique ID
  std::string name_;                    //!< Name of material
  vector<int> nuclide_;                 //!< Indices in nuclides vector
  vector<int> element_;                 //!< Indices in elements vector
  NCrystalMat ncrystal_mat_;            //!< NCrystal material object
  xt::xtensor<double, 1> atom_density_; //!< Nuclide atom density in [atom/b-cm]
  double density_;                      //!< Total atom density in [atom/b-cm]
  double density_gpcc_;                 //!< Total atom density in [g/cm^3]
  double volume_ {-1.0};                //!< Volume in [cm^3]
  vector<bool> p0_; //!< Indicate which nuclides are to be treated with
                    //!< iso-in-lab scattering

  // To improve performance of tallying, we store an array (direct address
  // table) that indicates for each nuclide in data::nuclides the index of the
  // corresponding nuclide in the nuclide_ vector. If it is not present in the
  // material, the entry is set to -1.
  vector<int> mat_nuclide_index_;

  // Thermal scattering tables
  vector<ThermalTable> thermal_tables_;

  unique_ptr<Bremsstrahlung> ttb_;

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
  int64_t index_;

  bool depletable_ {false}; //!< Is the material depletable?
  bool fissionable_ {
    false}; //!< Does this material contain fissionable nuclides
  //! \brief Default temperature for cells containing this material.
  //!
  //! A negative value indicates no default temperature was specified.
  double temperature_ {-1};
};

//==============================================================================
class Stochastic_Media {
public:
  //----------------------------------------------------------------------------
  // Constructors, destructors, factory functions

  explicit Stochastic_Media(pugi::xml_node cell_node);
  Stochastic_Media() {};
  virtual ~Stochastic_Media();

  // Accessors
  //! Get name
  //! \return Material name
  const std::string& name() const { return name_; }

  //! Set name
  void set_name(const std::string& name) { name_ = name; }

  //! Get ID of material
  //! \return ID of material
  int32_t id() const { return id_; }

  //! Assign a unique ID to the material
  //! \param[in] Unique ID to assign. A value of -1 indicates that an ID
  //!   should be automatically assigned.
  void set_id(int32_t id_);

  //! Get the particle radius
  double radius() const { return radius_; }

  //! Get the packing fraction of the particle material
  double pf() const { return pf_; }

  //! Get the material of the particle
  vector<int32_t> particle_mat() const { return particle_mat_; }

  //! Get the material of the matrix
  vector<int32_t> matrix_mat() const { return matrix_mat_; }

  //----------------------------------------------------------------------------
  // Data
  int32_t id_ {C_NONE}; //!< Unique ID
  std::string name_;    //!< Name of stochastic_media
  // Parameters: the  particle radius, now this method only support sphere
  double radius_;
  // The packing fraction of the particle material
  double pf_;
  // The material of the particle
  vector<int32_t> particle_mat_;
  // The material of the matrix
  vector<int32_t> matrix_mat_;
  //----------------------------------------------------------------------------

  virtual double sample_chord_length(uint64_t* seed_ptr) = 0;

private:
  // Private data members
  int64_t index_;
};

class CLS_Media : public Stochastic_Media {
public:
  explicit CLS_Media(pugi::xml_node cell_node);
  CLS_Media() {};

  double sample_chord_length(uint64_t* seed_ptr) override
  {
    double matrix_mean_chord = 4 / 3 * radius_ * (1 - pf_) / pf_;
    return -matrix_mean_chord * std::log(prn(seed_ptr));
  }
  double sample_particle_length(uint64_t* seed_ptr)
  {
    double cos_value = sqrt(prn(seed_ptr));
    return 2 * radius_ * cos_value;
  }
};

//==============================================================================
// Non-member functions
//==============================================================================

//! Calculate Sternheimer adjustment factor
double sternheimer_adjustment(const vector<double>& f,
  const vector<double>& e_b_sq, double e_p_sq, double n_conduction,
  double log_I, double tol, int max_iter);

//! Calculate density effect correction
double density_effect(const vector<double>& f, const vector<double>& e_b_sq,
  double e_p_sq, double n_conduction, double rho, double E, double tol,
  int max_iter);

//! Read material data from materials.xml
void read_materials_xml();

//! Read material data XML node
//! \param[in] root node of materials XML element
void read_materials_xml(pugi::xml_node root);

void free_memory_material();

} // namespace openmc
#endif // OPENMC_MATERIAL_H
