#ifndef OPENMC_PHOTON_H
#define OPENMC_PHOTON_H

#include "openmc/endf.h"
#include "openmc/memory.h" // for unique_ptr
#include "openmc/particle.h"
#include "openmc/vector.h"

#include <gsl/gsl>
#include <hdf5.h>
#include "xtensor/xtensor.hpp"

#include <string>
#include <unordered_map>
#include <utility> // for pair

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
  xsfloat n_electrons;
  xsfloat binding_energy;
  xt::xtensor<xsfloat, 1> cross_section;

  // Transition data
  int n_transitions;
  xt::xtensor<int, 2> transition_subshells;
  xt::xtensor<xsfloat, 1> transition_energy;
  xt::xtensor<xsfloat, 1> transition_probability;
};

class PhotonInteraction {
public:
  // Constructors/destructor
  PhotonInteraction(hid_t group);
  ~PhotonInteraction();

  // Methods
  void calculate_xs(Particle& p) const;

  void compton_scatter(xsfloat alpha, bool doppler, xsfloat* alpha_out,
    xsfloat* mu, int* i_shell, uint64_t* seed) const;

  xsfloat rayleigh_scatter(xsfloat alpha, uint64_t* seed) const;

  void pair_production(xsfloat alpha, xsfloat* E_electron, xsfloat* E_positron,
    xsfloat* mu_electron, xsfloat* mu_positron, uint64_t* seed) const;

  void atomic_relaxation(const ElectronSubshell& shell, Particle& p) const;

  // Data members
  std::string name_; //!< Name of element, e.g. "Zr"
  int Z_; //!< Atomic number
  gsl::index index_; //!< Index in global elements vector

  // Microscopic cross sections
  xt::xtensor<xsfloat, 1> energy_;
  xt::xtensor<xsfloat, 1> coherent_;
  xt::xtensor<xsfloat, 1> incoherent_;
  xt::xtensor<xsfloat, 1> photoelectric_total_;
  xt::xtensor<xsfloat, 1> pair_production_total_;
  xt::xtensor<xsfloat, 1> pair_production_electron_;
  xt::xtensor<xsfloat, 1> pair_production_nuclear_;
  xt::xtensor<xsfloat, 1> heating_;

  // Form factors
  Tabulated1D incoherent_form_factor_;
  Tabulated1D coherent_int_form_factor_;
  Tabulated1D coherent_anomalous_real_;
  Tabulated1D coherent_anomalous_imag_;

  // Photoionization and atomic relaxation data
  std::unordered_map<int, int> shell_map_; //!< Given a shell designator, e.g. 3, this
                                           //!< dictionary gives an index in shells_
  vector<ElectronSubshell> shells_;

  // Compton profile data
  xt::xtensor<xsfloat, 2> profile_pdf_;
  xt::xtensor<xsfloat, 2> profile_cdf_;
  xt::xtensor<xsfloat, 1> binding_energy_;
  xt::xtensor<xsfloat, 1> electron_pdf_;

  // Stopping power data
  xsfloat I_; // mean excitation energy
  xt::xtensor<int, 1> n_electrons_;
  xt::xtensor<xsfloat, 1> ionization_energy_;
  xt::xtensor<xsfloat, 1> stopping_power_radiative_;

  // Bremsstrahlung scaled DCS
  xt::xtensor<xsfloat, 2> dcs_;

private:
  void compton_doppler(xsfloat alpha, xsfloat mu, xsfloat* E_out, int* i_shell,
                       uint64_t* seed) const;
};

//==============================================================================
// Non-member functions
//==============================================================================

std::pair<xsfloat, xsfloat> klein_nishina(xsfloat alpha, uint64_t* seed);

void free_memory_photon();

//==============================================================================
// Global variables
//==============================================================================

namespace data {

extern xt::xtensor<xsfloat, 1> compton_profile_pz; //! Compton profile momentum grid

//! Photon interaction data for each element
extern std::unordered_map<std::string, int> element_map;
extern vector<unique_ptr<PhotonInteraction>> elements;

} // namespace data

} // namespace openmc

#endif // OPENMC_PHOTON_H
