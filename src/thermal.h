#ifndef OPENMC_THERMAL_SCATTERING_H
#define OPENMC_THERMAL_SCATTERING_H

#include <cstddef>
#include <string>
#include <vector>

#include "xtensor/xtensor.hpp"

#include "hdf5_interface.h"
#include "nuclide.h"

namespace openmc {

// Secondary energy mode for S(a,b) inelastic scattering
// TODO: Convert to enum
constexpr int SAB_SECONDARY_EQUAL  {0}; // Equally-likely outgoing energy bins
constexpr int SAB_SECONDARY_SKEWED {1}; // Skewed outgoing energy bins
constexpr int SAB_SECONDARY_CONT   {2}; // Continuous, linear-linear interpolation

// Elastic mode for S(a,b) elastic scattering
// TODO: Convert to enum
constexpr int SAB_ELASTIC_INCOHERENT {3}; // Incoherent elastic scattering
constexpr int SAB_ELASTIC_COHERENT   {4}; // Coherent elastic scattering (Bragg edges)

class ThermalData {
public:
  ThermalData(hid_t group, int secondary_mode);

  // Sample an outgoing energy and angle
  void sample(const NuclideMicroXS* micro_xs, double E_in,
              double* E_out, double* mu);
private:
  struct DistEnergySab {
    std::size_t n_e_out;
    xt::xtensor<double, 1> e_out;
    xt::xtensor<double, 1> e_out_pdf;
    xt::xtensor<double, 1> e_out_cdf;
    xt::xtensor<double, 2> mu;
  };

  // Threshold for thermal scattering treatment (usually ~4 eV)
  double threshold_inelastic_;
  double threshold_elastic_ {0.0};

  // Inelastic scattering data
  int inelastic_mode_;            // secondary mode (equal/skewed/continuous)
  std::size_t n_inelastic_e_in_;  // # of incoming E for inelastic
  std::size_t n_inelastic_e_out_; // # of outgoing E for inelastic
  std::size_t n_inelastic_mu_;    // # of outgoing angles for inelastic
  std::vector<double> inelastic_e_in_;
  std::vector<double> inelastic_sigma_;

  // The following are used only if secondary_mode is 0 or 1
  xt::xtensor<double, 2> inelastic_e_out_;
  xt::xtensor<double, 3> inelastic_mu_;

  // The following is used only if secondary_mode is 3
  // The different implementation is necessary because the continuous
  // representation has a variable number of outgoing energy points for each
  // incoming energy
  std::vector<DistEnergySab> inelastic_data_; // One for each Ein

  // Elastic scattering data
  int elastic_mode_;   // elastic mode (discrete/exact)
  std::size_t n_elastic_e_in_; // # of incoming E for elastic
  std::size_t n_elastic_mu_;   // # of outgoing angles for elastic
  std::vector<double> elastic_e_in_;
  std::vector<double> elastic_P_;
  xt::xtensor<double, 2> elastic_mu_;

  friend class ThermalScattering;
};

class ThermalScattering {
public:
  ThermalScattering(hid_t group, const std::vector<double>& temperature, int method,
                    double tolerance, const double* minmax);

  void calculate_xs(double E, double sqrtkT, int* i_temp, double* elastic,
                    double* inelastic) const ;
  bool has_nuclide(const char* name) const;

  std::string name_;    // name of table, e.g. "c_H_in_H2O"
  double awr_;      // weight of nucleus in neutron masses
  std::vector<double> kTs_;  // temperatures in eV (k*T)
  std::vector<std::string> nuclides_; // List of valid nuclides

  // cross sections and distributions at each temperature
  std::vector<ThermalData> data_;
};

} // namespace openmc

#endif // OPENMC_THERMAL_SCATTERING_H
