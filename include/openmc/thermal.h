#ifndef OPENMC_THERMAL_H
#define OPENMC_THERMAL_H

#include <cstddef>
#include <string>
#include <vector>

#include "xtensor/xtensor.hpp"

#include "openmc/hdf5_interface.h"
#include "openmc/nuclide.h"

namespace openmc {

//==============================================================================
// Constants
//==============================================================================

// Secondary energy mode for S(a,b) inelastic scattering
// TODO: Convert to enum
constexpr int SAB_SECONDARY_EQUAL  {0}; // Equally-likely outgoing energy bins
constexpr int SAB_SECONDARY_SKEWED {1}; // Skewed outgoing energy bins
constexpr int SAB_SECONDARY_CONT   {2}; // Continuous, linear-linear interpolation

// Elastic mode for S(a,b) elastic scattering
// TODO: Convert to enum
constexpr int SAB_ELASTIC_INCOHERENT {3}; // Incoherent elastic scattering
constexpr int SAB_ELASTIC_COHERENT   {4}; // Coherent elastic scattering (Bragg edges)

//==============================================================================
//! Secondary angle-energy data for thermal neutron scattering at a single
//! temperature
//==============================================================================

class ThermalData {
public:
  ThermalData(hid_t group, int secondary_mode);

  // Sample an outgoing energy and angle
  void sample(const NuclideMicroXS* micro_xs, double E_in,
              double* E_out, double* mu);
private:
  //! Secondary energy/angle distributions for inelastic thermal scattering
  //! collisions which utilize a continuous secondary energy representation.
  struct DistEnergySab {
    std::size_t n_e_out; //!< Number of outgoing energies
    xt::xtensor<double, 1> e_out;     //!< Outgoing energies
    xt::xtensor<double, 1> e_out_pdf; //!< Probability density function
    xt::xtensor<double, 1> e_out_cdf; //!< Cumulative distribution function
    xt::xtensor<double, 2> mu; //!< Equiprobable angles at each outgoing energy
  };

  //! Upper threshold for incoherent inelastic scattering (usually ~4 eV)
  double threshold_inelastic_;
  //! Upper threshold for coherent/incoherent elastic scattering
  double threshold_elastic_ {0.0};

  // Inelastic scattering data
  int inelastic_mode_;            //!< distribution type (equal/skewed/continuous)
  std::size_t n_inelastic_e_in_;  //!< number of incoming E for inelastic
  std::size_t n_inelastic_e_out_; //!< number of outgoing E for inelastic
  std::size_t n_inelastic_mu_;    //!< number of outgoing angles for inelastic
  std::vector<double> inelastic_e_in_; //!< incoming E grid for inelastic
  std::vector<double> inelastic_sigma_; //!< inelastic scattering cross section

  // The following are used only for equal/skewed distributions
  xt::xtensor<double, 2> inelastic_e_out_;
  xt::xtensor<double, 3> inelastic_mu_;

  // The following is used only for continuous S(a,b) distributions. The
  // different implementation is necessary because the continuous representation
  // has a variable number of outgoing energy points for each incoming energy
  std::vector<DistEnergySab> inelastic_data_; //!< Secondary angle-energy at
                                              //!< each incoming energy

  // Elastic scattering data
  int elastic_mode_;   //!< type of elastic (incoherent/coherent)
  std::size_t n_elastic_e_in_; //!< number of incoming E for elastic
  std::size_t n_elastic_mu_;   //!< number of outgoing angles for elastic
  std::vector<double> elastic_e_in_; //!< incoming E grid for elastic
  std::vector<double> elastic_P_; //!< elastic scattering cross section
  xt::xtensor<double, 2> elastic_mu_; //!< equi-probable angles at each incoming E

  // ThermalScattering needs access to private data members
  friend class ThermalScattering;
};

//==============================================================================
//! Data for thermal neutron scattering, typically off light isotopes in
//! moderating materials such as water, graphite, BeO, etc.
//==============================================================================

class ThermalScattering {
public:
  ThermalScattering(hid_t group, const std::vector<double>& temperature, int method,
                    double tolerance, const double* minmax);

  //! Determine inelastic/elastic cross section at given energy
  //!
  //! \param[in] E incoming energy in [eV]
  //! \param[in] sqrtkT square-root of temperature multipled by Boltzmann's constant
  //! \param[out] i_temp corresponding temperature index
  //! \param[out] elastic Thermal elastic scattering cross section
  //! \param[out] inelastic Thermal inelastic scattering cross section
  void calculate_xs(double E, double sqrtkT, int* i_temp, double* elastic,
                    double* inelastic) const;

  //! Determine whether table applies to a particular nuclide
  //!
  //! \param[in] name Name of the nuclide, e.g., "H1"
  //! \return Whether table applies to the nuclide
  bool has_nuclide(const char* name) const;

  // Sample an outgoing energy and angle
  void sample(const NuclideMicroXS* micro_xs, double E_in,
              double* E_out, double* mu);

  double threshold() const { return data_[0].threshold_inelastic_; }

  std::string name_; //!< name of table, e.g. "c_H_in_H2O"
  double awr_;       //!< weight of nucleus in neutron masses
  std::vector<double> kTs_;  //!< temperatures in [eV] (k*T)
  std::vector<std::string> nuclides_; //!< Valid nuclides

  //! cross sections and distributions at each temperature
  std::vector<ThermalData> data_;
};

//==============================================================================
// Fortran compatibility functions
//==============================================================================

extern "C" {
  ThermalScattering* sab_from_hdf5(hid_t group, const double* temperature,
    int n, int method, double tolerance, const double* minmax);
  void sab_calculate_xs(ThermalScattering* data, double E, double sqrtkT,
    int* i_temp, double* elastic, double* inelastic);
  void sab_free(ThermalScattering* data);
  bool sab_has_nuclide(ThermalScattering* data, const char* name);
  void sab_sample(ThermalScattering* data, const NuclideMicroXS* micro_xs,
    double E_in, double* E_out, double* mu);
  double sab_threshold(ThermalScattering* data);
}

} // namespace openmc

#endif // OPENMC_THERMAL_H
