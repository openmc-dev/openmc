#ifndef OPENMC_THERMAL_H
#define OPENMC_THERMAL_H

#include <cstddef>
#include <unordered_map>

#include "xtensor/xtensor.hpp"

#include "openmc/angle_energy.h"
#include "openmc/endf.h"
#include "openmc/hdf5_interface.h"
#include "openmc/memory.h"
#include "openmc/particle.h"
#include "openmc/string.h"
#include "openmc/vector.h"

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

class ThermalScattering;

// TODO put thermal scattering on the GPU

namespace data {
extern std::unordered_map<std::string, int> thermal_scatt_map;
extern vector<unique_ptr<ThermalScattering>> thermal_scatt;
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
  void calculate_xs(xsfloat E, xsfloat* elastic, xsfloat* inelastic) const;

  //! Sample an outgoing energy and angle
  //
  //! \param[in] micro_xs Microscopic cross sections
  //! \param[in] E_in Incident neutron energy in [eV]
  //! \param[out] E_out Outgoing neutron energy in [eV]
  //! \param[out] mu Outgoing scattering angle cosine
  //! \param[inout] seed Pseudorandom seed pointer
  void sample(const NuclideMicroXS& micro_xs, xsfloat E_in,
              xsfloat* E_out, xsfloat* mu, uint64_t* seed);
private:
  struct Reaction {
    // Default constructor
    Reaction() { }

    // Data members
    unique_ptr<Function1D> xs; //!< Cross section
    unique_ptr<AngleEnergy>
      distribution; //!< Secondary angle-energy distribution
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
  ThermalScattering(hid_t group, const vector<xsfloat>& temperature);

  //! Determine inelastic/elastic cross section at given energy
  //!
  //! \param[in] E incoming energy in [eV]
  //! \param[in] sqrtkT square-root of temperature multipled by Boltzmann's constant
  //! \param[out] i_temp corresponding temperature index
  //! \param[out] elastic Thermal elastic scattering cross section
  //! \param[out] inelastic Thermal inelastic scattering cross section
  //! \param[inout] seed Pseudorandom seed pointer
  void calculate_xs(xsfloat E, xsfloat sqrtkT, int* i_temp, xsfloat* elastic,
                    xsfloat* inelastic, uint64_t* seed) const;

  //! Determine whether table applies to a particular nuclide
  //!
  //! \param[in] name Name of the nuclide, e.g., "H1"
  //! \return Whether table applies to the nuclide
  bool has_nuclide(const char* name) const;

  // Sample an outgoing energy and angle
  void sample(const NuclideMicroXS& micro_xs, xsfloat E_in,
              xsfloat* E_out, xsfloat* mu);

  std::string name_; //!< name of table, e.g. "c_H_in_H2O"
  double awr_;       //!< weight of nucleus in neutron masses
  double energy_max_; //!< maximum energy for thermal scattering in [eV]
  vector<xsfloat> kTs_;           //!< temperatures in [eV] (k*T)
  vector<string> nuclides_; //!< Valid nuclides

  //! cross sections and distributions at each temperature
  vector<ThermalData> data_;
};

void free_memory_thermal();

} // namespace openmc

#endif // OPENMC_THERMAL_H
