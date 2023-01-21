#ifndef OPENMC_PHOTON_H
#define OPENMC_PHOTON_H

#include "openmc/endf.h"
#include "openmc/particle.h"
#include "openmc/serialize.h"
#include "openmc/vector.h"

#include <gsl/gsl>
#include <hdf5.h>
#include "xtensor/xtensor.hpp"

#include <memory> // for unique_ptr
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
  struct Transition {
    int primary_subshell;
    int secondary_subshell;
    double energy;
    double probability;
  };

  // Constructors
  ElectronSubshell() { };

  int index_subshell; //!< index in SUBSHELLS
  int threshold;
  double n_electrons;
  double binding_energy;
  gsl::span<const double> cross_section;
  gsl::span<const Transition> transitions;
};

class PhotonInteraction {
public:
  // Constructors/destructor
  PhotonInteraction(hid_t group);
  ~PhotonInteraction();

  // Methods
  #pragma omp declare target
  void calculate_xs(Particle& p) const;

  void compton_scatter(double alpha, bool doppler, double* alpha_out,
    double* mu, int* i_shell, uint64_t* seed) const;

  double rayleigh_scatter(double alpha, uint64_t* seed) const;

  void pair_production(double alpha, double* E_electron, double* E_positron,
    double* mu_electron, double* mu_positron, uint64_t* seed) const;

  void atomic_relaxation(int i_shell, Particle& p) const;
  #pragma omp end declare target

  void copy_to_device();
  void release_from_device();

  // Data members
  std::string name_; //!< Name of element, e.g. "Zr"
  int Z_; //!< Atomic number
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

  double* device_energy_ {nullptr};
  double* device_coherent_ {nullptr};
  double* device_incoherent_ {nullptr};
  double* device_pair_production_total_ {nullptr};

  // Form factors
  Tabulated1D incoherent_form_factor_;
  Tabulated1D coherent_int_form_factor_;
  Tabulated1D coherent_anomalous_real_;
  Tabulated1D coherent_anomalous_imag_;

  Tabulated1DFlat incoherent_form_factor() const;
  Tabulated1DFlat coherent_int_form_factor() const;

  // Photoionization and atomic relaxation data
  vector<ElectronSubshell> shells_;

  // Compton profile data
  xt::xtensor<double, 2> profile_pdf_;
  xt::xtensor<double, 2> profile_cdf_;
  xt::xtensor<double, 1> binding_energy_;
  xt::xtensor<double, 1> electron_pdf_;

  size_t n_profile_;
  double* device_profile_pdf_ {nullptr};
  double* device_profile_cdf_ {nullptr};
  double* device_binding_energy_ {nullptr};
  double* device_electron_pdf_ {nullptr};
  double profile_pdf(gsl::index i, gsl::index j) const;
  double profile_cdf(gsl::index i, gsl::index j) const;

  // Stopping power data
  double I_; // mean excitation energy
  xt::xtensor<int, 1> n_electrons_;
  xt::xtensor<double, 1> ionization_energy_;
  xt::xtensor<double, 1> stopping_power_radiative_;

  // Bremsstrahlung scaled DCS
  xt::xtensor<double, 2> dcs_;

  // Constant data
  static constexpr int MAX_STACK_SIZE =
    7; //!< maximum possible size of atomic relaxation stack

private:
  void compton_doppler(double alpha, double mu, double* E_out, int* i_shell,
                       uint64_t* seed) const;

  int calc_max_stack_size() const;
  int calc_helper(std::unordered_map<int, int>& visited, int i_shell) const;

  // Backend storage
  vector<double> cross_section_;
  vector<ElectronSubshell::Transition> transitions_;

  DataBuffer buffer_; // form factors
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

#pragma omp declare target
extern double* compton_profile_pz; //! Compton profile momentum grid
extern size_t compton_profile_pz_size;
#pragma omp end declare target

//! Photon interaction data for each element
extern std::unordered_map<std::string, int> element_map;
#pragma omp declare target
extern PhotonInteraction* elements;
extern size_t elements_size;
#pragma omp end declare target
extern size_t elements_capacity;

} // namespace data

} // namespace openmc

#endif // OPENMC_PHOTON_H
