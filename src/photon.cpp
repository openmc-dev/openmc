#include "openmc/photon.h"

#include "openmc/array.h"
#include "openmc/bremsstrahlung.h"
#include "openmc/constants.h"
#include "openmc/distribution_multi.h"
#include "openmc/hdf5_interface.h"
#include "openmc/message_passing.h"
#include "openmc/nuclide.h"
#include "openmc/particle.h"
#include "openmc/physics.h"
#include "openmc/random_dist.h"
#include "openmc/random_lcg.h"
#include "openmc/search.h"
#include "openmc/settings.h"

#include "openmc/tensor.h"

#include <cmath>
#include <fmt/core.h>
#include <limits>
#include <stdexcept>
#include <tuple> // for tie

namespace openmc {

constexpr int PhotonInteraction::MAX_STACK_SIZE;

//==============================================================================
// Global variables
//==============================================================================

namespace data {

tensor::Tensor<double> compton_profile_pz;

std::unordered_map<std::string, int> element_map;
vector<unique_ptr<PhotonInteraction>> elements;

} // namespace data

//==============================================================================
// PhotonInteraction implementation
//==============================================================================

PhotonInteraction::PhotonInteraction(hid_t group)
{
  // Set index of element in global vector
  index_ = data::elements.size();

  // Get name of nuclide from group, removing leading '/'
  name_ = object_name(group).substr(1);
  data::element_map[name_] = index_;

  // Get atomic number
  read_attribute(group, "Z", Z_);

  // Determine number of energies and read energy grid
  read_dataset(group, "energy", energy_);

  // Read coherent scattering
  hid_t rgroup = open_group(group, "coherent");
  read_dataset(rgroup, "xs", coherent_);

  hid_t dset = open_dataset(rgroup, "integrated_scattering_factor");
  coherent_int_form_factor_ = Tabulated1D {dset};
  close_dataset(dset);

  if (object_exists(group, "anomalous_real")) {
    dset = open_dataset(rgroup, "anomalous_real");
    coherent_anomalous_real_ = Tabulated1D {dset};
    close_dataset(dset);
  }

  if (object_exists(group, "anomalous_imag")) {
    dset = open_dataset(rgroup, "anomalous_imag");
    coherent_anomalous_imag_ = Tabulated1D {dset};
    close_dataset(dset);
  }
  close_group(rgroup);

  // Read incoherent scattering
  rgroup = open_group(group, "incoherent");
  read_dataset(rgroup, "xs", incoherent_);
  dset = open_dataset(rgroup, "scattering_factor");
  incoherent_form_factor_ = Tabulated1D {dset};
  close_dataset(dset);
  close_group(rgroup);

  // Read pair production
  if (object_exists(group, "pair_production_electron")) {
    rgroup = open_group(group, "pair_production_electron");
    read_dataset(rgroup, "xs", pair_production_electron_);
    close_group(rgroup);
  } else {
    pair_production_electron_ = tensor::zeros_like(energy_);
  }

  // Read pair production
  if (object_exists(group, "pair_production_nuclear")) {
    rgroup = open_group(group, "pair_production_nuclear");
    read_dataset(rgroup, "xs", pair_production_nuclear_);
    close_group(rgroup);
  } else {
    pair_production_nuclear_ = tensor::zeros_like(energy_);
  }

  // Read photoelectric
  rgroup = open_group(group, "photoelectric");
  read_dataset(rgroup, "xs", photoelectric_total_);
  close_group(rgroup);

  // Read heating
  if (object_exists(group, "heating")) {
    rgroup = open_group(group, "heating");
    read_dataset(rgroup, "xs", heating_);
    close_group(rgroup);
  } else {
    heating_ = tensor::zeros_like(energy_);
  }

  // Read subshell photoionization cross section and atomic relaxation data
  rgroup = open_group(group, "subshells");
  vector<std::string> designators;
  read_attribute(rgroup, "designators", designators);
  auto n_shell = designators.size();
  if (n_shell == 0) {
    throw std::runtime_error {
      "Photoatomic data for " + name_ + " does not have subshell data."};
  }

  shells_.resize(n_shell);
  cross_sections_ = tensor::zeros<double>({energy_.size(), n_shell});

  // Create mapping from designator to index
  std::unordered_map<int, int> shell_map;
  for (int i = 0; i < n_shell; ++i) {
    const auto& designator {designators[i]};

    int j = 1;
    for (const auto& subshell : SUBSHELLS) {
      if (designator == subshell) {
        shell_map[j] = i;
        shells_[i].index_subshell = j;
        break;
      }
      ++j;
    }
  }
  shell_map[0] = -1;

  for (int i = 0; i < n_shell; ++i) {
    const auto& designator {designators[i]};
    auto& shell {shells_[i]};

    // TODO: Move to ElectronSubshell constructor

    hid_t tgroup = open_group(rgroup, designator.c_str());

    // Read binding energy if atomic relaxation data is present
    if (attribute_exists(tgroup, "binding_energy")) {
      has_atomic_relaxation_ = true;
      read_attribute(tgroup, "binding_energy", shell.binding_energy);
    }

    // Read subshell cross section
    tensor::Tensor<double> xs;
    dset = open_dataset(tgroup, "xs");
    read_attribute(dset, "threshold_idx", shell.threshold);
    close_dataset(dset);
    read_dataset(tgroup, "xs", xs);

    auto cross_section =
      cross_sections_.slice(tensor::range(static_cast<size_t>(shell.threshold),
                              cross_sections_.shape(0)),
        i);
    cross_section = tensor::where(xs > 0, tensor::log(xs), 0);

    if (settings::atomic_relaxation && object_exists(tgroup, "transitions")) {
      // Determine dimensions of transitions
      dset = open_dataset(tgroup, "transitions");
      auto dims = object_shape(dset);
      close_dataset(dset);

      int n_transition = dims[0];
      if (n_transition > 0) {
        tensor::Tensor<double> matrix;
        read_dataset(tgroup, "transitions", matrix);

        // Transition probability normalization
        double norm =
          tensor::Tensor<double>(matrix.slice(tensor::all, 3)).sum();

        shell.transitions.resize(n_transition);
        for (int j = 0; j < n_transition; ++j) {
          auto& transition = shell.transitions[j];
          transition.primary_subshell = shell_map.at(matrix(j, 0));
          transition.secondary_subshell = shell_map.at(matrix(j, 1));
          transition.energy = matrix(j, 2);
          transition.probability = matrix(j, 3) / norm;
        }
      }
    }
    close_group(tgroup);
  }
  close_group(rgroup);

  // Check the maximum size of the atomic relaxation stack
  auto max_size = this->calc_max_stack_size();
  if (max_size > MAX_STACK_SIZE && mpi::master) {
    warning(fmt::format("The subshell vacancy stack in atomic relaxation can "
                        "grow up to {}, but the stack size limit is set to {}.",
      max_size, MAX_STACK_SIZE));
  }

  // Determine number of electron shells
  rgroup = open_group(group, "compton_profiles");

  // Read electron shell PDF and binding energies
  read_dataset(rgroup, "num_electrons", electron_pdf_);
  electron_pdf_ /= electron_pdf_.sum();
  read_dataset(rgroup, "binding_energy", binding_energy_);

  // Read Compton profiles
  read_dataset(rgroup, "J", profile_pdf_);

  // Get Compton profile momentum grid
  if (data::compton_profile_pz.size() == 0) {
    read_dataset(rgroup, "pz", data::compton_profile_pz);
  }
  close_group(rgroup);

  // Map Compton subshell data to atomic relaxation data by finding the
  // subshell with the equivalent binding energy
  if (settings::atomic_relaxation && has_atomic_relaxation_) {
    auto is_close = [](double a, double b) {
      return std::abs(a - b) / a < FP_REL_PRECISION;
    };
    subshell_map_ = tensor::Tensor<int>(binding_energy_.shape(), -1);
    for (int i = 0; i < binding_energy_.size(); ++i) {
      double E_b = binding_energy_[i];
      if (i < n_shell && is_close(E_b, shells_[i].binding_energy)) {
        subshell_map_[i] = i;
      } else {
        for (int j = 0; j < n_shell; ++j) {
          if (is_close(E_b, shells_[j].binding_energy)) {
            subshell_map_[i] = j;
            break;
          }
        }
      }
    }
  }

  // Create Compton profile CDF
  auto n_profile = data::compton_profile_pz.size();
  auto n_shell_compton = profile_pdf_.shape(0);
  if (n_profile < 2) {
    throw std::runtime_error {
      "At least two points are required in a Compton profile."};
  }
  profile_cdf_ = tensor::Tensor<double>({n_shell_compton, n_profile});
  profile_tail_slope_ = tensor::Tensor<double>({n_shell_compton});
  profile_negative_mass_ = tensor::Tensor<double>({n_shell_compton});
  if (n_shell_compton > SUBSHELLS.size()) {
    throw std::runtime_error {"Photoatomic data for element " + name_ +
                              " has more Compton profiles than supported "
                              "electron subshells."};
  }
  for (int i = 0; i < n_shell_compton; ++i) {
    double c = 0.0;
    profile_cdf_(i, 0) = 0.0;
    for (int j = 0; j < n_profile - 1; ++j) {
      c += 0.5 *
           (data::compton_profile_pz(j + 1) - data::compton_profile_pz(j)) *
           (profile_pdf_(i, j) + profile_pdf_(i, j + 1));
      profile_cdf_(i, j + 1) = c;
    }

    // Extrapolate the profile beyond the tabulated grid linearly on a
    // log-linear scale. The normalization includes the extrapolated tail.
    double pz_last = data::compton_profile_pz(n_profile - 1);
    double pz_prev = data::compton_profile_pz(n_profile - 2);
    double profile_last = profile_pdf_(i, n_profile - 1);
    double profile_prev = profile_pdf_(i, n_profile - 2);
    if (!(pz_last > pz_prev) || !(profile_last > 0.0) ||
        !(profile_prev > 0.0)) {
      throw std::runtime_error {"The final two points of the Compton profile "
                                "for element " +
                                name_ + " are not valid for extrapolation."};
    }
    double slope = std::log(profile_last / profile_prev) / (pz_last - pz_prev);
    if (!std::isfinite(slope) || slope >= 0.0) {
      throw std::runtime_error {"The final two values of the Compton profile "
                                "for element " +
                                name_ + " do not form a decreasing tail."};
    }
    profile_tail_slope_(i) = slope;
    double norm = 2.0 * (c - profile_last / slope);
    if (!std::isfinite(norm) || norm <= 0.0) {
      throw std::runtime_error {"The Compton profile for element " + name_ +
                                " has an invalid normalization."};
    }
    for (int j = 0; j < n_profile; ++j) {
      profile_pdf_(i, j) /= norm;
      profile_cdf_(i, j) /= norm;
    }
  }
  for (int i = 0; i < n_shell_compton; ++i) {
    profile_negative_mass_(i) = this->compton_profile_cdf(i, FINE_STRUCTURE);
  }

  // Calculate total pair production
  pair_production_total_ = pair_production_nuclear_ + pair_production_electron_;

  if (settings::electron_treatment == ElectronTreatment::TTB) {
    // Read bremsstrahlung scaled DCS
    rgroup = open_group(group, "bremsstrahlung");
    read_dataset(rgroup, "dcs", dcs_);
    auto n_e = dcs_.shape(0);
    auto n_k = dcs_.shape(1);

    // Get energy grids used for bremsstrahlung DCS and for stopping powers
    tensor::Tensor<double> electron_energy;
    read_dataset(rgroup, "electron_energy", electron_energy);
    if (data::ttb_k_grid.size() == 0) {
      read_dataset(rgroup, "photon_energy", data::ttb_k_grid);
    }

    // Get data used for density effect correction
    read_dataset(rgroup, "num_electrons", n_electrons_);
    read_dataset(rgroup, "ionization_energy", ionization_energy_);
    read_attribute(rgroup, "I", I_);
    close_group(rgroup);

    // Truncate the bremsstrahlung data at the cutoff energy
    int photon = ParticleType::photon().transport_index();
    const auto& E {electron_energy};
    double cutoff = settings::energy_cutoff[photon];
    if (cutoff > E(0)) {
      size_t i_grid = lower_bound_index(
        E.cbegin(), E.cend(), settings::energy_cutoff[photon]);

      // calculate interpolation factor
      double f = (std::log(cutoff) - std::log(E(i_grid))) /
                 (std::log(E(i_grid + 1)) - std::log(E(i_grid)));

      // Interpolate bremsstrahlung DCS at the cutoff energy and truncate
      tensor::Tensor<double> dcs({n_e - i_grid, n_k});
      for (int i = 0; i < n_k; ++i) {
        double y = std::exp(
          std::log(dcs_(i_grid, i)) +
          f * (std::log(dcs_(i_grid + 1, i)) - std::log(dcs_(i_grid, i))));
        tensor::View<double> col_i = dcs.slice(tensor::all, i);
        col_i(0) = y;
        for (int j = i_grid + 1; j < n_e; ++j) {
          col_i(j - i_grid) = dcs_(j, i);
        }
      }
      dcs_ = dcs;

      tensor::Tensor<double> frst({static_cast<size_t>(1)});
      frst(0) = cutoff;
      tensor::Tensor<double> rest(electron_energy.slice(
        tensor::range(i_grid + 1, electron_energy.size())));
      electron_energy = tensor::concatenate(frst, rest);
    }

    // Set incident particle energy grid
    if (data::ttb_e_grid.size() == 0) {
      data::ttb_e_grid = electron_energy;
    }

    // Calculate the radiative stopping power
    stopping_power_radiative_ =
      tensor::Tensor<double>({data::ttb_e_grid.size()});
    for (int i = 0; i < data::ttb_e_grid.size(); ++i) {
      // Integrate over reduced photon energy
      double c = 0.0;
      for (int j = 0; j < data::ttb_k_grid.size() - 1; ++j) {
        c += 0.5 * (dcs_(i, j + 1) + dcs_(i, j)) *
             (data::ttb_k_grid(j + 1) - data::ttb_k_grid(j));
      }
      double e = data::ttb_e_grid(i);

      // Square of the ratio of the speed of light to the velocity of the
      // charged particle
      double beta_sq = e * (e + 2.0 * MASS_ELECTRON_EV) /
                       ((e + MASS_ELECTRON_EV) * (e + MASS_ELECTRON_EV));

      stopping_power_radiative_(i) = Z_ * Z_ / beta_sq * e * c;
    }
  }

  // Take logarithm of energies and cross sections since they are log-log
  // interpolated. Note that cross section libraries converted from ACE files
  // represent zero as exp(-500) to avoid log-log interpolation errors. For
  // values below exp(-499) we store the log as -900, for which exp(-900)
  // evaluates to zero.
  double limit = std::exp(-499.0);
  energy_ = tensor::log(energy_);
  coherent_ = tensor::where(coherent_ > limit, tensor::log(coherent_), -900.0);
  incoherent_ =
    tensor::where(incoherent_ > limit, tensor::log(incoherent_), -900.0);
  photoelectric_total_ = tensor::where(
    photoelectric_total_ > limit, tensor::log(photoelectric_total_), -900.0);
  pair_production_total_ = tensor::where(pair_production_total_ > limit,
    tensor::log(pair_production_total_), -900.0);
  heating_ = tensor::where(heating_ > limit, tensor::log(heating_), -900.0);
}

PhotonInteraction::~PhotonInteraction()
{
  data::element_map.erase(name_);
}

int PhotonInteraction::calc_max_stack_size() const
{
  // Table to store solutions to sub-problems
  std::unordered_map<int, int> visited;

  // Find the maximum possible size of the stack used to store holes created
  // during atomic relaxation, checking over every subshell the initial hole
  // could be in
  int max_size = 0;
  for (int i_shell = 0; i_shell < shells_.size(); ++i_shell) {
    max_size = std::max(max_size, this->calc_helper(visited, i_shell));
  }
  return max_size;
}

int PhotonInteraction::calc_helper(
  std::unordered_map<int, int>& visited, int i_shell) const
{
  // No transitions for this subshell, so this is the only shell in the stack
  const auto& shell {shells_[i_shell]};
  if (shell.transitions.empty()) {
    return 1;
  }

  // Check the table to see if the maximum stack size has already been
  // calculated for this shell
  auto it = visited.find(i_shell);
  if (it != visited.end()) {
    return it->second;
  }

  int max_size = 0;
  for (const auto& transition : shell.transitions) {
    // If this is a non-radiative transition two vacancies are created and
    // the stack grows by one; if this is a radiative transition only one
    // vacancy is created and the stack size stays the same
    int size = 0;
    if (transition.secondary_subshell != -1) {
      size = this->calc_helper(visited, transition.secondary_subshell) + 1;
    }
    size =
      std::max(size, this->calc_helper(visited, transition.primary_subshell));
    max_size = std::max(max_size, size);
  }
  visited[i_shell] = max_size;
  return max_size;
}

void PhotonInteraction::compton_scatter(double alpha, bool doppler,
  double* alpha_out, double* mu, int* i_shell, uint64_t* seed) const
{
  double form_factor_xmax = 0.0;
  while (true) {
    // Sample Klein-Nishina distribution for trial energy and angle
    std::tie(*alpha_out, *mu) = klein_nishina(alpha, seed);

    // Note that the parameter used here does not correspond exactly to the
    // momentum transfer q in ENDF-102 Eq. (27.2). Rather, this is the
    // parameter as defined by Hubbell, where the actual data comes from
    double x =
      MASS_ELECTRON_EV / PLANCK_C * alpha * std::sqrt(0.5 * (1.0 - *mu));

    // Calculate S(x, Z) and S(x_max, Z)
    double form_factor_x = incoherent_form_factor_(x);
    if (form_factor_xmax == 0.0) {
      form_factor_xmax =
        incoherent_form_factor_(MASS_ELECTRON_EV / PLANCK_C * alpha);
    }

    // Perform rejection on form factor
    if (prn(seed) < form_factor_x / form_factor_xmax) {
      if (doppler) {
        double E_out;
        this->compton_doppler(alpha, *mu, &E_out, i_shell, seed);
        *alpha_out = E_out / MASS_ELECTRON_EV;
      } else {
        *i_shell = -1;
      }
      break;
    }
  }
}

double PhotonInteraction::compton_profile_cdf(int i_shell, double pz) const
{
  if (pz <= 0.0)
    return 0.0;

  auto n = data::compton_profile_pz.size();
  double pz_last = data::compton_profile_pz(n - 1);
  double c;
  if (pz >= pz_last) {
    c = profile_cdf_(i_shell, n - 1) + detail::compton_profile_tail_integral(pz,
                                         pz_last, profile_pdf_(i_shell, n - 1),
                                         profile_tail_slope_(i_shell));
  } else {
    int i = lower_bound_index(
      data::compton_profile_pz.cbegin(), data::compton_profile_pz.cend(), pz);
    double pz_l = data::compton_profile_pz(i);
    double pz_r = data::compton_profile_pz(i + 1);
    double p_l = profile_pdf_(i_shell, i);
    double p_r = profile_pdf_(i_shell, i + 1);
    double c_l = profile_cdf_(i_shell, i);
    double slope = (p_r - p_l) / (pz_r - pz_l);
    double delta = pz - pz_l;
    c = c_l + p_l * delta + 0.5 * slope * delta * delta;
  }
  return std::min(0.5, c);
}

double PhotonInteraction::invert_compton_profile_cdf(
  int i_shell, double c) const
{
  auto n = data::compton_profile_pz.size();
  double integral = c;
  double c_last = profile_cdf_(i_shell, n - 1);
  if (integral >= c_last) {
    // Invert the log-linear extrapolated tail using Kaltiaisenaho Eq. (3.123).
    return detail::invert_compton_profile_tail(integral - c_last,
      data::compton_profile_pz(n - 1), profile_pdf_(i_shell, n - 1),
      profile_tail_slope_(i_shell));
  }

  // Invert the piecewise-linear tabulated profile (Kaltiaisenaho Eq. 3.126).
  // The rationalized quadratic root used below is equivalent to that equation
  // but remains well-conditioned when the profile slope is small.
  tensor::View<const double> cdf_shell = profile_cdf_.slice(i_shell);
  int i = lower_bound_index(cdf_shell.cbegin(), cdf_shell.cend(), integral);
  double pz_l = data::compton_profile_pz(i);
  double pz_r = data::compton_profile_pz(i + 1);
  double p_l = profile_pdf_(i_shell, i);
  double p_r = profile_pdf_(i_shell, i + 1);
  double c_l = profile_cdf_(i_shell, i);
  if (p_l == p_r) {
    return pz_l + (integral - c_l) / p_l;
  }

  double slope = (p_r - p_l) / (pz_r - pz_l);
  double delta_c = integral - c_l;
  double discriminant = p_l * p_l + 2.0 * slope * delta_c;
  double denominator = p_l + std::sqrt(std::max(0.0, discriminant));
  return pz_l + 2.0 * delta_c / denominator;
}

PhotonInteraction::ShellKinematics PhotonInteraction::compton_shell_kinematics(
  double alpha, double mu, double E, int i_shell) const
{
  ShellKinematics kinematics {};
  double E_b = binding_energy_(i_shell);
  if (E <= E_b)
    return kinematics;

  // Kaltiaisenaho Eq. (3.73): substitute E' = E - E_b in the RIA
  // kinematic relation to obtain the upper bound on the allowed pz interval.
  kinematics.pz_max = -FINE_STRUCTURE * (E_b - (E - E_b) * alpha * (1.0 - mu)) /
                      std::sqrt(2.0 * E * (E - E_b) * (1.0 - mu) + E_b * E_b);
  if (kinematics.pz_max <= -FINE_STRUCTURE)
    return kinematics;

  // Eq. (3.118), using profile_negative_mass_ = K_i(1/alpha): the
  // kinematically accessible mass is the integral from -1/alpha to pz_max.
  double c_negative = profile_negative_mass_(i_shell);
  kinematics.c_limit =
    this->compton_profile_cdf(i_shell, std::abs(kinematics.pz_max));
  kinematics.profile_mass =
    c_negative + std::copysign(kinematics.c_limit, kinematics.pz_max);
  return kinematics;
}

bool PhotonInteraction::sample_compton_momentum(double alpha, double mu,
  double E, int i_shell, const ShellKinematics& kinematics, double* E_out,
  uint64_t* seed) const
{
  double c_negative = profile_negative_mass_(i_shell);
  double pz;
  // Inverse-transform sampling of the signed pz distribution, following
  // Kaltiaisenaho Eqs. (3.120)-(3.126). The tabulated profile is symmetric,
  // so its negative branch is obtained by reflecting the half-profile CDF.
  if (kinematics.pz_max < 0.0) {
    double c = kinematics.c_limit + prn(seed) * kinematics.profile_mass;
    pz = -this->invert_compton_profile_cdf(i_shell, c);
  } else {
    double c = prn(seed) * kinematics.profile_mass;
    if (c < c_negative) {
      pz = -this->invert_compton_profile_cdf(i_shell, c_negative - c);
    } else {
      pz = this->invert_compton_profile_cdf(i_shell, c - c_negative);
    }
  }

  double energy_ratio = detail::compton_energy_ratio(alpha, mu, pz);
  double max_energy_ratio = 1.0 - binding_energy_(i_shell) / E;
  if (!std::isfinite(energy_ratio) || energy_ratio <= 0.0)
    return false;

  double energy_tolerance = 16.0 * std::numeric_limits<double>::epsilon() *
                            std::max(1.0, max_energy_ratio);
  if (energy_ratio > max_energy_ratio + energy_tolerance)
    return false;

  energy_ratio = std::min(energy_ratio, max_energy_ratio);
  *E_out = energy_ratio * E;

  // Kaltiaisenaho Eq. (3.127): account for the E'/E factor in the
  // approximate RIA DDCS after solving the scattered-photon energy.
  return prn(seed) <= energy_ratio;
}

bool PhotonInteraction::compton_doppler_conditional(double alpha, double mu,
  double E, double* E_out, int* i_shell, uint64_t* seed) const
{
  array<ShellKinematics, SUBSHELLS.size()> shell_data;
  array<double, SUBSHELLS.size()> shell_cdf;
  double shell_pmf_norm = 0.0;

  // Form the shell PMF in Eq. (3.116), f_i times the accessible profile mass.
  // This is algebraically equivalent to repeated shell rejection in Eq.
  // (3.119).
  for (int i = 0; i < electron_pdf_.size(); ++i) {
    shell_data[i] = this->compton_shell_kinematics(alpha, mu, E, i);
    shell_pmf_norm += electron_pdf_(i) * shell_data[i].profile_mass;
    shell_cdf[i] = shell_pmf_norm;
  }
  if (shell_pmf_norm == 0.0)
    return false;

  // The conditional shell PMF avoids repeated rejection when the accessible
  // profile mass is small. Retain a bound to protect against degenerate
  // momentum/energy sampling and roundoff.
  constexpr int MAX_SAMPLES = 100000;
  for (int attempt = 0; attempt < MAX_SAMPLES; ++attempt) {
    double rn = prn(seed) * shell_pmf_norm;
    int shell;
    for (shell = 0; shell < electron_pdf_.size(); ++shell) {
      if (rn < shell_cdf[shell])
        break;
    }
    *i_shell = shell;

    if (this->sample_compton_momentum(
          alpha, mu, E, shell, shell_data[shell], E_out, seed))
      return true;
  }
  return false;
}

void PhotonInteraction::compton_doppler(
  double alpha, double mu, double* E_out, int* i_shell, uint64_t* seed) const
{
  // Implements the approximate RIA Doppler-broadening algorithm in Sec. 3.4.8
  // of T. Kaltiaisenaho, "Implementing a photon physics model in Serpent 2"
  // (2016), https://aaltodoc.aalto.fi/handle/123456789/21004.
  // First use Kaltiaisenaho's shell-rejection procedure (Eqs. 3.116-3.119),
  // which usually accepts quickly. If it does not, sample its equivalent
  // conditional shell PMF to bound work for near-forward scattering.
  constexpr int N_FAST_SAMPLES = 2;

  double E = alpha * MASS_ELECTRON_EV;
  int shell = 0;
  for (int attempt = 0; attempt < N_FAST_SAMPLES; ++attempt) {
    // Propose shell i according to occupancy f_i (first step of Eq. 3.119).
    double rn = prn(seed);
    double c = 0.0;
    for (shell = 0; shell < electron_pdf_.size(); ++shell) {
      c += electron_pdf_(shell);
      if (rn < c)
        break;
    }

    auto kinematics = this->compton_shell_kinematics(alpha, mu, E, shell);
    if (kinematics.profile_mass <= 0.0)
      continue;

    // Accept with the accessible Compton-profile mass (Eq. 3.119).
    if (prn(seed) >= kinematics.profile_mass)
      continue;

    if (this->sample_compton_momentum(
          alpha, mu, E, shell, kinematics, E_out, seed)) {
      *i_shell = shell;
      return;
    }
  }

  *i_shell = shell;
  if (this->compton_doppler_conditional(alpha, mu, E, E_out, i_shell, seed))
    return;

  // No shell/momentum sample was accepted within the iteration budget.
  // Fall back to the free-electron Compton energy for the last sampled
  // shell rather than looping indefinitely.
  *E_out = alpha / (1.0 + alpha * (1.0 - mu)) * MASS_ELECTRON_EV;
}

void PhotonInteraction::calculate_xs(Particle& p) const
{
  // Perform binary search on the element energy grid in order to determine
  // which points to interpolate between
  int n_grid = energy_.size();
  double log_E = std::log(p.E());
  int i_grid;
  if (log_E <= energy_[0]) {
    i_grid = 0;
  } else if (log_E > energy_(n_grid - 1)) {
    i_grid = n_grid - 2;
  } else {
    // We use upper_bound_index here because sometimes photons are created with
    // energies that exactly match a grid point
    i_grid = upper_bound_index(energy_.cbegin(), energy_.cend(), log_E);
  }

  // check for case where two energy points are the same
  if (energy_(i_grid) == energy_(i_grid + 1))
    ++i_grid;

  // calculate interpolation factor
  double f =
    (log_E - energy_(i_grid)) / (energy_(i_grid + 1) - energy_(i_grid));

  auto& xs {p.photon_xs(index_)};
  xs.index_grid = i_grid;
  xs.interp_factor = f;

  // Calculate microscopic coherent cross section
  xs.coherent = std::exp(
    coherent_(i_grid) + f * (coherent_(i_grid + 1) - coherent_(i_grid)));

  // Calculate microscopic incoherent cross section
  xs.incoherent = std::exp(
    incoherent_(i_grid) + f * (incoherent_(i_grid + 1) - incoherent_(i_grid)));

  // Calculate microscopic photoelectric cross section
  xs.photoelectric = 0.0;
  tensor::View<const double> xs_lower = cross_sections_.slice(i_grid);
  tensor::View<const double> xs_upper = cross_sections_.slice(i_grid + 1);

  for (int i = 0; i < xs_upper.size(); ++i)
    if (xs_lower(i) != 0)
      xs.photoelectric +=
        std::exp(xs_lower(i) + f * (xs_upper(i) - xs_lower(i)));

  // Calculate microscopic pair production cross section
  xs.pair_production = std::exp(
    pair_production_total_(i_grid) +
    f * (pair_production_total_(i_grid + 1) - pair_production_total_(i_grid)));

  // Calculate microscopic total cross section
  xs.total =
    xs.coherent + xs.incoherent + xs.photoelectric + xs.pair_production;
  xs.last_E = p.E();
}

double PhotonInteraction::rayleigh_scatter(double alpha, uint64_t* seed) const
{
  double mu;
  while (true) {
    // Determine maximum value of x^2
    double x2_max = std::pow(MASS_ELECTRON_EV / PLANCK_C * alpha, 2);

    // Determine F(x^2_max, Z)
    double F_max = coherent_int_form_factor_(x2_max);

    // Sample cumulative distribution
    double F = prn(seed) * F_max;

    // Determine x^2 corresponding to F
    const auto& x {coherent_int_form_factor_.x()};
    const auto& y {coherent_int_form_factor_.y()};
    int i = lower_bound_index(y.cbegin(), y.cend(), F);
    double r = (F - y[i]) / (y[i + 1] - y[i]);
    double x2 = x[i] + r * (x[i + 1] - x[i]);

    // Calculate mu
    mu = 1.0 - 2.0 * x2 / x2_max;

    if (prn(seed) < 0.5 * (1.0 + mu * mu))
      break;
  }
  return mu;
}

void PhotonInteraction::pair_production(double alpha, double* E_electron,
  double* E_positron, double* mu_electron, double* mu_positron,
  uint64_t* seed) const
{
  constexpr double r[] {122.81, 73.167, 69.228, 67.301, 64.696, 61.228, 57.524,
    54.033, 50.787, 47.851, 46.373, 45.401, 44.503, 43.815, 43.074, 42.321,
    41.586, 40.953, 40.524, 40.256, 39.756, 39.144, 38.462, 37.778, 37.174,
    36.663, 35.986, 35.317, 34.688, 34.197, 33.786, 33.422, 33.068, 32.740,
    32.438, 32.143, 31.884, 31.622, 31.438, 31.142, 30.950, 30.758, 30.561,
    30.285, 30.097, 29.832, 29.581, 29.411, 29.247, 29.085, 28.930, 28.721,
    28.580, 28.442, 28.312, 28.139, 27.973, 27.819, 27.675, 27.496, 27.285,
    27.093, 26.911, 26.705, 26.516, 26.304, 26.108, 25.929, 25.730, 25.577,
    25.403, 25.245, 25.100, 24.941, 24.790, 24.655, 24.506, 24.391, 24.262,
    24.145, 24.039, 23.922, 23.813, 23.712, 23.621, 23.523, 23.430, 23.331,
    23.238, 23.139, 23.048, 22.967, 22.833, 22.694, 22.624, 22.545, 22.446,
    22.358, 22.264};

  // The reduced screening radius r is the ratio of the screening radius to
  // the Compton wavelength of the electron, where the screening radius is
  // obtained under the assumption that the Coulomb field of the nucleus is
  // exponentially screened by atomic electrons. This allows us to use a
  // simplified atomic form factor and analytical approximations of the
  // screening functions in the pair production DCS instead of computing the
  // screening functions numerically. The reduced screening radii above for
  // Z = 1-99 come from F. Salvat, J. M. Fernández-Varea, and J. Sempau,
  // "PENELOPE-2011: A Code System for Monte Carlo Simulation of Electron and
  // Photon Transport," OECD-NEA, Issy-les-Moulineaux, France (2011).

  // Compute the high-energy Coulomb correction
  double a = Z_ / FINE_STRUCTURE;
  double c =
    a * a *
    (1.0 / (1.0 + a * a) + 0.202059 +
      a * a *
        (-0.03693 +
          a * a *
            (0.00835 +
              a * a *
                (-0.00201 +
                  a * a * (0.00049 + a * a * (-0.00012 + a * a * 0.00003))))));

  // The analytical approximation of the DCS underestimates the cross section
  // at low energies. The correction factor f compensates for this.
  double q = std::sqrt(2.0 / alpha);
  double f = q * (-0.1774 - 12.10 * a + 11.18 * a * a) +
             q * q * (8.523 + 73.26 * a - 44.41 * a * a) +
             q * q * q * (-13.52 - 121.1 * a + 96.41 * a * a) +
             q * q * q * q * (8.946 + 62.05 * a - 63.41 * a * a);

  // Calculate phi_1(1/2) and phi_2(1/2). The unnormalized PDF for the reduced
  // energy is given by p = 2*(1/2 - e)^2*phi_1(e) + phi_2(e), where phi_1 and
  // phi_2 are non-negative and maximum at e = 1/2.
  double b = 2.0 * r[Z_] / alpha;
  double t1 = 2.0 * std::log(1.0 + b * b);
  double t2 = b * std::atan(1.0 / b);
  double t3 = b * b * (4.0 - 4.0 * t2 - 3.0 * std::log(1.0 + 1.0 / (b * b)));
  double t4 = 4.0 * std::log(r[Z_]) - 4.0 * c + f;
  double phi1_max = 7.0 / 3.0 - t1 - 6.0 * t2 - t3 + t4;
  double phi2_max = 11.0 / 6.0 - t1 - 3.0 * t2 + 0.5 * t3 + t4;

  // To aid sampling, the unnormalized PDF can be expressed as
  // p = u_1*U_1(e)*pi_1(e) + u_2*U_2(e)*pi_2(e), where pi_1 and pi_2 are
  // normalized PDFs on the interval (e_min, e_max) from which values of e can
  // be sampled using the inverse transform method, and
  // U_1 = phi_1(e)/phi_1(1/2) and U_2 = phi_2(e)/phi_2(1/2) are valid
  // rejection functions. The reduced energy can now be sampled using a
  // combination of the composition and rejection methods.
  double u1 = 2.0 / 3.0 * std::pow(0.5 - 1.0 / alpha, 2) * phi1_max;
  double u2 = phi2_max;
  double e;
  while (true) {
    double rn = prn(seed);

    // Sample the index i in (1, 2) using the point probabilities
    // p(1) = u_1/(u_1 + u_2) and p(2) = u_2/(u_1 + u_2)
    int i;
    if (prn(seed) < u1 / (u1 + u2)) {
      i = 1;

      // Sample e from pi_1 using the inverse transform method
      e = rn >= 0.5
            ? 0.5 + (0.5 - 1.0 / alpha) * std::pow(2.0 * rn - 1.0, 1.0 / 3.0)
            : 0.5 - (0.5 - 1.0 / alpha) * std::pow(1.0 - 2.0 * rn, 1.0 / 3.0);
    } else {
      i = 2;

      // Sample e from pi_2 using the inverse transform method
      e = 1.0 / alpha + (0.5 - 1.0 / alpha) * 2.0 * rn;
    }

    // Calculate phi_i(e) and deliver e if rn <= U_i(e)
    b = r[Z_] / (2.0 * alpha * e * (1.0 - e));
    t1 = 2.0 * std::log(1.0 + b * b);
    t2 = b * std::atan(1.0 / b);
    t3 = b * b * (4.0 - 4.0 * t2 - 3.0 * std::log(1.0 + 1.0 / (b * b)));
    if (i == 1) {
      double phi1 = 7.0 / 3.0 - t1 - 6.0 * t2 - t3 + t4;
      if (prn(seed) <= phi1 / phi1_max)
        break;
    } else {
      double phi2 = 11.0 / 6.0 - t1 - 3.0 * t2 + 0.5 * t3 + t4;
      if (prn(seed) <= phi2 / phi2_max)
        break;
    }
  }

  // Compute the kinetic energy of the electron and the positron
  *E_electron = (alpha * e - 1.0) * MASS_ELECTRON_EV;
  *E_positron = (alpha * (1.0 - e) - 1.0) * MASS_ELECTRON_EV;

  // Sample the scattering angle of the electron. The cosine of the polar
  // angle of the direction relative to the incident photon is sampled from
  // p(mu) = C/(1 - beta*mu)^2 using the inverse transform method.
  double beta =
    std::sqrt(*E_electron * (*E_electron + 2.0 * MASS_ELECTRON_EV)) /
    (*E_electron + MASS_ELECTRON_EV);
  double rn = uniform_distribution(-1., 1., seed);
  *mu_electron = (rn + beta) / (rn * beta + 1.0);

  // Sample the scattering angle of the positron
  beta = std::sqrt(*E_positron * (*E_positron + 2.0 * MASS_ELECTRON_EV)) /
         (*E_positron + MASS_ELECTRON_EV);
  rn = uniform_distribution(-1., 1., seed);
  *mu_positron = (rn + beta) / (rn * beta + 1.0);
}

void PhotonInteraction::atomic_relaxation(int i_shell, Particle& p) const
{
  // Return if no atomic relaxation data is present or if the binding energy is
  // larger than the incident particle energy
  if (!has_atomic_relaxation_ || shells_[i_shell].binding_energy > p.E())
    return;

  // Stack for unprocessed holes left by transitioning electrons
  int n_holes = 0;
  array<int, MAX_STACK_SIZE> holes;

  // Push the initial hole onto the stack
  holes[n_holes++] = i_shell;

  while (n_holes > 0) {
    // Pop the next hole off the stack
    int i_hole = holes[--n_holes];
    const auto& shell {shells_[i_hole]};

    // If no transitions, assume fluorescent photon from captured free electron
    if (shell.transitions.empty()) {
      Direction u = isotropic_direction(p.current_seed());
      double E = shell.binding_energy;
      p.create_secondary(p.wgt(), u, E, ParticleType::photon());
      continue;
    }

    // Sample transition
    double c = -prn(p.current_seed());
    int i_trans;
    for (i_trans = 0; i_trans < shell.transitions.size(); ++i_trans) {
      c += shell.transitions[i_trans].probability;
      if (c > 0)
        break;
    }
    const auto& transition = shell.transitions[i_trans];

    // Sample angle isotropically
    Direction u = isotropic_direction(p.current_seed());

    // Push the hole created by the electron transitioning to the photoelectron
    // hole onto the stack
    holes[n_holes++] = transition.primary_subshell;

    if (transition.secondary_subshell != -1) {
      // Non-radiative transition -- Auger/Coster-Kronig effect

      // Push the hole left by emitted auger electron onto the stack
      holes[n_holes++] = transition.secondary_subshell;

      // Process Auger electron at the photon collision site.
      process_charged_secondary(
        p, u, transition.energy, ParticleType::electron());
    } else {
      // Radiative transition -- get X-ray energy

      // Create fluorescent photon
      p.create_secondary(p.wgt(), u, transition.energy, ParticleType::photon());
    }
  }
}

//==============================================================================
// Non-member functions
//==============================================================================

double detail::compton_profile_tail_integral(
  double pz, double pz_last, double profile_last, double slope)
{
  return profile_last * std::expm1(slope * (pz - pz_last)) / slope;
}

double detail::invert_compton_profile_tail(
  double integral, double pz_last, double profile_last, double slope)
{
  return pz_last + std::log1p(slope * integral / profile_last) / slope;
}

double detail::compton_energy_ratio(double alpha, double mu, double pz)
{
  if (pz == 0.0)
    return 1.0 / (1.0 + alpha * (1.0 - mu));

  double momentum = pz / FINE_STRUCTURE;
  double momentum_sq = momentum * momentum;
  double f = 1.0 + alpha * (1.0 - mu);
  double a = momentum_sq - f * f;
  double b = 2.0 * (f - momentum_sq * mu);
  double c = momentum_sq - 1.0;
  double discriminant = b * b - 4.0 * a * c;
  double discriminant_tolerance = 16.0 *
                                  std::numeric_limits<double>::epsilon() *
                                  (b * b + std::abs(4.0 * a * c));
  if (discriminant < -discriminant_tolerance)
    return std::numeric_limits<double>::quiet_NaN();
  discriminant = std::max(0.0, discriminant);

  double root1;
  double root2;
  if (std::abs(a) < 1.0e-14 * (std::abs(b) + std::abs(c))) {
    if (b == 0.0)
      return std::numeric_limits<double>::quiet_NaN();
    root1 = -c / b;
    root2 = root1;
  } else {
    double sqrt_discriminant = std::sqrt(discriminant);
    double q = -0.5 * (b + std::copysign(sqrt_discriminant, b));
    root1 = q / a;
    root2 = q == 0.0 ? (-b + sqrt_discriminant) / (2.0 * a) : c / q;
  }

  double root_min = std::numeric_limits<double>::infinity();
  double root_max = -std::numeric_limits<double>::infinity();
  if (std::isfinite(root1) && root1 > 0.0) {
    root_min = root1;
    root_max = root1;
  }
  if (std::isfinite(root2) && root2 > 0.0) {
    root_min = std::min(root_min, root2);
    root_max = std::max(root_max, root2);
  }
  if (!std::isfinite(root_min))
    return std::numeric_limits<double>::quiet_NaN();

  double energy_ratio = pz < 0.0 ? root_min : root_max;
  double free_electron_ratio = 1.0 / f;
  double tolerance = 16.0 * std::numeric_limits<double>::epsilon() *
                     std::max(1.0, free_electron_ratio);
  if ((pz < 0.0 && energy_ratio > free_electron_ratio + tolerance) ||
      (pz > 0.0 && energy_ratio < free_electron_ratio - tolerance)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return energy_ratio;
}

std::pair<double, double> klein_nishina(double alpha, uint64_t* seed)
{
  double alpha_out, mu;
  double beta = 1.0 + 2.0 * alpha;
  if (alpha < 3.0) {
    // Kahn's rejection method
    double t = beta / (beta + 8.0);
    double x;
    while (true) {
      if (prn(seed) < t) {
        // Left branch of flow chart
        double r = uniform_distribution(0.0, 2.0, seed);
        x = 1.0 + alpha * r;
        if (prn(seed) < 4.0 / x * (1.0 - 1.0 / x)) {
          mu = 1 - r;
          break;
        }
      } else {
        // Right branch of flow chart
        x = beta / (1.0 + 2.0 * alpha * prn(seed));
        mu = 1.0 + (1.0 - x) / alpha;
        if (prn(seed) < 0.5 * (mu * mu + 1.0 / x))
          break;
      }
    }
    alpha_out = alpha / x;

  } else {
    // Koblinger's direct method
    double gamma = 1.0 - std::pow(beta, -2);
    double s =
      prn(seed) * (4.0 / alpha + 0.5 * gamma +
                    (1.0 - (1.0 + beta) / (alpha * alpha)) * std::log(beta));
    if (s <= 2.0 / alpha) {
      // For first term, x = 1 + 2ar
      // Therefore, a' = a/(1 + 2ar)
      alpha_out = alpha / (1.0 + 2.0 * alpha * prn(seed));
    } else if (s <= 4.0 / alpha) {
      // For third term, x = beta/(1 + 2ar)
      // Therefore, a' = a(1 + 2ar)/beta
      alpha_out = alpha * (1.0 + 2.0 * alpha * prn(seed)) / beta;
    } else if (s <= 4.0 / alpha + 0.5 * gamma) {
      // For fourth term, x = 1/sqrt(1 - gamma*r)
      // Therefore, a' = a*sqrt(1 - gamma*r)
      alpha_out = alpha * std::sqrt(1.0 - gamma * prn(seed));
    } else {
      // For third term, x = beta^r
      // Therefore, a' = a/beta^r
      alpha_out = alpha / std::pow(beta, prn(seed));
    }

    // Calculate cosine of scattering angle based on basic relation
    mu = 1.0 + 1.0 / alpha - 1.0 / alpha_out;
  }
  return {alpha_out, mu};
}

void free_memory_photon()
{
  data::elements.clear();
  data::compton_profile_pz.resize({0});
  data::ttb_e_grid.resize({0});
  data::ttb_k_grid.resize({0});
}

} // namespace openmc
