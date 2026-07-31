#ifndef OPENMC_PHOTON_H
#define OPENMC_PHOTON_H

#include "openmc/endf.h"
#include "openmc/memory.h" // for unique_ptr
#include "openmc/particle.h"
#include "openmc/vector.h"

#include "openmc/tensor.h"
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
  int64_t index_;    //!< Index in global elements vector

  // Microscopic cross sections
  tensor::Tensor<double> energy_;
  tensor::Tensor<double> coherent_;
  tensor::Tensor<double> incoherent_;
  tensor::Tensor<double> photoelectric_total_;
  tensor::Tensor<double> pair_production_total_;
  tensor::Tensor<double> pair_production_electron_;
  tensor::Tensor<double> pair_production_nuclear_;
  tensor::Tensor<double> heating_;

  // Form factors
  Tabulated1D incoherent_form_factor_;
  Tabulated1D coherent_int_form_factor_;
  Tabulated1D coherent_anomalous_real_;
  Tabulated1D coherent_anomalous_imag_;

  // Photoionization and atomic relaxation data. Subshell cross sections are
  // stored separately to improve memory access pattern when calculating the
  // total cross section
  vector<ElectronSubshell> shells_;
  tensor::Tensor<double> cross_sections_;

  // Compton profile data
  tensor::Tensor<double> profile_pdf_;
  tensor::Tensor<double> profile_cdf_;
  tensor::Tensor<double> profile_tail_slope_;
  tensor::Tensor<double> profile_negative_mass_; //!< Mass from -1/alpha to 0
  tensor::Tensor<double> binding_energy_;
  tensor::Tensor<double> electron_pdf_;

  // Map subshells from Compton profile data obtained from Biggs et al,
  // "Hartree-Fock Compton profiles for the elements" to ENDF/B atomic
  // relaxation data
  tensor::Tensor<int> subshell_map_;

  // Stopping power data
  double I_; // mean excitation energy
  tensor::Tensor<int> n_electrons_;
  tensor::Tensor<double> ionization_energy_;
  tensor::Tensor<double> stopping_power_radiative_;

  // Bremsstrahlung scaled DCS
  tensor::Tensor<double> dcs_;

  // Whether atomic relaxation data is present
  bool has_atomic_relaxation_ {false};

  // Constant data
  static constexpr int MAX_STACK_SIZE =
    7; //!< maximum possible size of atomic relaxation stack
private:
  struct ShellKinematics {
    double pz_max;       //!< Upper bound in Kaltiaisenaho Eq. (3.73)
    double c_limit;      //!< Half-profile integral K_i(|pz_max|), Eq. (3.117)
    double profile_mass; //!< Accessible profile mass, Eq. (3.118)
  };

  void compton_doppler(
    double alpha, double mu, double* E_out, int* i_shell, uint64_t* seed) const;

  //! Determine pz_max and the accessible profile mass (Eqs. 3.73, 3.118)
  ShellKinematics compton_shell_kinematics(
    double alpha, double mu, double E, int i_shell) const;

  //! Sample pz and E' for a selected shell (Eqs. 3.120-3.127)
  bool sample_compton_momentum(double alpha, double mu, double E, int i_shell,
    const ShellKinematics& kinematics, double* E_out, uint64_t* seed) const;

  //! Sample from the shell PMF in Kaltiaisenaho Eq. (3.116)
  bool compton_doppler_conditional(double alpha, double mu, double E,
    double* E_out, int* i_shell, uint64_t* seed) const;

  //! Evaluate K_i(pz), the normalized half-profile integral (Eq. 3.117)
  double compton_profile_cdf(int i_shell, double pz) const;

  //! Invert K_i using the inverse transforms of Eqs. (3.123) and (3.126)
  double invert_compton_profile_cdf(int i_shell, double c) const;

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

namespace detail {

//! Integrate an exponentially extrapolated Compton-profile tail
double compton_profile_tail_integral(
  double pz, double pz_last, double profile_last, double slope);

//! Invert an exponentially extrapolated Compton-profile tail integral
double invert_compton_profile_tail(
  double integral, double pz_last, double profile_last, double slope);

//! Calculate the outgoing-to-incident energy ratio for signed electron momentum
double compton_energy_ratio(double alpha, double mu, double pz);

} // namespace detail

std::pair<double, double> klein_nishina(double alpha, uint64_t* seed);

void free_memory_photon();

//==============================================================================
// Global variables
//==============================================================================

namespace data {

extern tensor::Tensor<double>
  compton_profile_pz; //! Compton profile momentum grid

//! Photon interaction data for each element
extern std::unordered_map<std::string, int> element_map;
extern vector<unique_ptr<PhotonInteraction>> elements;

} // namespace data

} // namespace openmc

#endif // OPENMC_PHOTON_H
