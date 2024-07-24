#ifndef OPENMC_PHOTON_H
#define OPENMC_PHOTON_H

#include "openmc/endf.h"
#include "openmc/memory.h" // for unique_ptr
#include "openmc/particle.h"
#include "openmc/vector.h"

#include "xtensor/xtensor.hpp"
#include <gsl/gsl-lite.hpp>
#include <hdf5.h>

#include <string>
#include <unordered_map>
#include <utility> // for pair

namespace openmc {

//==============================================================================
//! Photon interaction data for a single element
//==============================================================================

class ElectronSubshell {
public:
  struct Transition {
    int primary_subshell;   //!< Index in shells_ of originating subshell
    int secondary_subshell; //!< Index in shells_ of Auger electron subshell
    double energy;          //!< Energy of transition
    double probability;     //!< Probability of transition between subshells
  };

  // Constructors
  ElectronSubshell() {};

  int index_subshell; //!< index in SUBSHELLS
  int threshold;
  double n_electrons;
  double binding_energy;
  vector<Transition> transitions;
};

class PhotonInteraction {
public:
  // Constructors/destructor
  PhotonInteraction(hid_t group);
  ~PhotonInteraction();

  // Methods
  void calculate_xs(Particle& p) const;

  void compton_scatter(double alpha, bool doppler, double* alpha_out,
    double* mu, int* i_shell, uint64_t* seed) const;

  double rayleigh_scatter(double alpha, uint64_t* seed) const;

  void pair_production(double alpha, double* E_electron, double* E_positron,
    double* mu_electron, double* mu_positron, uint64_t* seed) const;

  void atomic_relaxation(int i_shell, Particle& p) const;

  // Data members
  std::string name_; //!< Name of element, e.g. "Zr"
  int Z_;            //!< Atomic number
  gsl::index index_; //!< Index in global elements vector

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

  // Photoionization and atomic relaxation data. Subshell cross sections are
  // stored separately to improve memory access pattern when calculating the
  // total cross section
  vector<ElectronSubshell> shells_;
  xt::xtensor<double, 2> cross_sections_;

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

  // Whether atomic relaxation data is present
  bool has_atomic_relaxation_ {false};

  // Constant data
  static constexpr int MAX_RELAXATION_STACK_SIZE =
    7; //!< maximum possible size of atomic relaxation stack
private:
  void compton_doppler(
    double alpha, double mu, double* E_out, int* i_shell, uint64_t* seed) const;

  //! Calculate the maximum size of the vacancy stack in atomic relaxation
  //
  //! These helper functions use the subshell transition data to calculate the
  //! maximum size the stack of unprocessed subshell vacancies can grow to for
  //! the given element while simulating the cascade of photons and electrons
  //! in atomic relaxation.
  int calc_max_stack_size() const;
  int calc_helper(std::unordered_map<int, int>& visited, int i_shell) const;
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

extern xt::xtensor<double, 1>
  compton_profile_pz; //! Compton profile momentum grid

//! Photon interaction data for each element
extern std::unordered_map<std::string, int> element_map;
extern vector<unique_ptr<PhotonInteraction>> elements;

} // namespace data

} // namespace openmc

#endif // OPENMC_PHOTON_H
