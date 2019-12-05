#ifndef OPENMC_PHOTON_H
#define OPENMC_PHOTON_H

#include "openmc/endf.h"
#include "openmc/particle.h"

#include <hdf5.h>
#include "xtensor/xtensor.hpp"

#include <string>
#include <unordered_map>
#include <utility> // for pair
#include <vector>

namespace openmc {

//==============================================================================
//! Photon interaction data for a single element
//==============================================================================

class ElectronSubshell {
public:
  // Constructors
  ElectronSubshell() { };

  int index_subshell;  //!< index in SUBSHELLS
  int threshold;
  double n_electrons;
  double binding_energy;
  xt::xtensor<double, 1> cross_section;

  // Transition data
  int n_transitions;
  xt::xtensor<int, 2> transition_subshells;
  xt::xtensor<double, 1> transition_energy;
  xt::xtensor<double, 1> transition_probability;
};

class PhotonInteraction {
public:
  // Constructors
  PhotonInteraction(hid_t group, int i_element);

  // Methods
  void calculate_xs(Particle& p) const;

  void compton_scatter(double alpha, bool doppler, double* alpha_out,
    double* mu, int* i_shell, uint64_t* seed) const;

  double rayleigh_scatter(double alpha, uint64_t* seed) const;

  void pair_production(double alpha, double* E_electron, double* E_positron,
    double* mu_electron, double* mu_positron, uint64_t* seed) const;

  void atomic_relaxation(const ElectronSubshell& shell, Particle& p) const;

  // Data members
  std::string name_; //!< Name of element, e.g. "Zr"
  int Z_; //!< Atomic number
  int i_element_; //!< Index in global elements vector

  // Microscopic cross sections
  xt::xtensor<double, 1> energy_;
  xt::xtensor<double, 1> coherent_;
  xt::xtensor<double, 1> incoherent_;
  xt::xtensor<double, 1> photoelectric_total_;
  xt::xtensor<double, 1> pair_production_total_;
  xt::xtensor<double, 1> pair_production_electron_;
  xt::xtensor<double, 1> pair_production_nuclear_;
  xt::xtensor<double, 1> heating_;

  // Form factors
  Tabulated1D incoherent_form_factor_;
  Tabulated1D coherent_int_form_factor_;
  Tabulated1D coherent_anomalous_real_;
  Tabulated1D coherent_anomalous_imag_;

  // Photoionization and atomic relaxation data
  std::unordered_map<int, int> shell_map_; //!< Given a shell designator, e.g. 3, this
                                           //!< dictionary gives an index in shells_
  std::vector<ElectronSubshell> shells_;

  // Compton profile data
  xt::xtensor<double, 2> profile_pdf_;
  xt::xtensor<double, 2> profile_cdf_;
  xt::xtensor<double, 1> binding_energy_;
  xt::xtensor<double, 1> electron_pdf_;

  // Stopping power data
  double I_; // mean excitation energy
  xt::xtensor<int, 1> n_electrons_;
  xt::xtensor<double, 1> ionization_energy_;
  xt::xtensor<double, 1> stopping_power_radiative_;

  // Bremsstrahlung scaled DCS
  xt::xtensor<double, 2> dcs_;

private:
  void compton_doppler(double alpha, double mu, double* E_out, int* i_shell,
                       uint64_t* seed) const;
};

//==============================================================================
// Non-member functions
//==============================================================================

std::pair<double, double> klein_nishina(double alpha, uint64_t* seed);

void free_memory_photon();

//==============================================================================
// Global variables
//==============================================================================

namespace data {

extern xt::xtensor<double, 1> compton_profile_pz; //! Compton profile momentum grid

//! Photon interaction data for each element
extern std::vector<PhotonInteraction> elements;
extern std::unordered_map<std::string, int> element_map;

} // namespace data

} // namespace openmc

#endif // OPENMC_PHOTON_H
