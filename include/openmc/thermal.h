#ifndef OPENMC_THERMAL_H
#define OPENMC_THERMAL_H

#include <cstddef>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include <iostream>
#include "xtensor/xtensor.hpp"

#include "openmc/angle_energy.h"
#include "openmc/endf.h"
#include "openmc/endf_flat.h"
#include "openmc/hdf5_interface.h"
#include "openmc/particle.h"
#include "openmc/secondary_flat.h"
#include "openmc/vector.h"

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

class ThermalScattering;

namespace data {
extern std::unordered_map<std::string, int> thermal_scatt_map;
extern std::vector<ThermalScattering> thermal_scatt;
#pragma omp declare target
extern ThermalScattering* device_thermal_scatt;
#pragma omp end declare target
}

//==============================================================================
//! Secondary angle-energy data for thermal neutron scattering at a single
//! temperature
//==============================================================================

class ThermalData {
public:
  ThermalData(hid_t group);

  //! Calculate the cross section
  //
  //! \param[in] E Incident neutron energy in [eV]
  //! \param[out] elastic Elastic scattering cross section in [b]
  //! \param[out] inelastic Inelastic scattering cross section in [b]
  #pragma omp declare target
  void calculate_xs(double E, double* elastic, double* inelastic) const;
  #pragma omp end declare target

  //! Sample an outgoing energy and angle
  //
  //! \param[in] micro_xs Microscopic cross sections
  //! \param[in] E_in Incident neutron energy in [eV]
  //! \param[out] E_out Outgoing neutron energy in [eV]
  //! \param[out] mu Outgoing scattering angle cosine
  //! \param[inout] seed Pseudorandom seed pointer
  #pragma omp declare target
  void sample(const NuclideMicroXS& micro_xs, double E_in,
              double* E_out, double* mu, uint64_t* seed);
  #pragma omp end declare target
private:
  struct Reaction {
    // Default constructor
    Reaction() { }

    // Data members
    std::unique_ptr<Function1DFlatContainer> xs; //!< Cross section
    std::unique_ptr<AngleEnergyFlatContainer> distribution; //!< Secondary angle-energy distribution
    Function1DFlatContainer* device_xs {nullptr};
    AngleEnergyFlatContainer* device_distribution {nullptr};
  };

  // Inelastic scattering data
  Reaction elastic_;
  Reaction inelastic_;

  // ThermalScattering needs access to private data members
  friend class ThermalScattering;
};

//==============================================================================
//! Data for thermal neutron scattering, typically off light isotopes in
//! moderating materials such as water, graphite, BeO, etc.
//==============================================================================

class ThermalScattering {
public:
  ThermalScattering(hid_t group, const std::vector<double>& temperature);

  //! Determine inelastic/elastic cross section at given energy
  //!
  //! \param[in] E incoming energy in [eV]
  //! \param[in] sqrtkT square-root of temperature multipled by Boltzmann's constant
  //! \param[out] i_temp corresponding temperature index
  //! \param[out] elastic Thermal elastic scattering cross section
  //! \param[out] inelastic Thermal inelastic scattering cross section
  //! \param[inout] seed Pseudorandom seed pointer
  #pragma omp declare target
  void calculate_xs(double E, double sqrtkT, int* i_temp, double* elastic,
                    double* inelastic, uint64_t* seed) const;
  #pragma omp end declare target

  //! Determine whether table applies to a particular nuclide
  //!
  //! \param[in] name Name of the nuclide, e.g., "H1"
  //! \return Whether table applies to the nuclide
  bool has_nuclide(const char* name) const;

  // Sample an outgoing energy and angle
  void sample(const NuclideMicroXS& micro_xs, double E_in,
              double* E_out, double* mu);

  void copy_to_device();
  void release_from_device();

  std::string name_; //!< name of table, e.g. "c_H_in_H2O"
  double awr_;       //!< weight of nucleus in neutron masses
  double energy_max_; //!< maximum energy for thermal scattering in [eV]
  vector<double> kTs_;  //!< temperatures in [eV] (k*T)
  std::vector<std::string> nuclides_; //!< Valid nuclides

  //! cross sections and distributions at each temperature
  vector<ThermalData> data_;
};

void free_memory_thermal();

} // namespace openmc

#endif // OPENMC_THERMAL_H
