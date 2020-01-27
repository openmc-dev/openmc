#include "openmc/photon.h"

#include "openmc/bremsstrahlung.h"
#include "openmc/constants.h"
#include "openmc/hdf5_interface.h"
#include "openmc/nuclide.h"
#include "openmc/particle.h"
#include "openmc/random_lcg.h"
#include "openmc/search.h"
#include "openmc/settings.h"

#include "xtensor/xbuilder.hpp"
#include "xtensor/xoperation.hpp"
#include "xtensor/xview.hpp"

#include <array>
#include <cmath>
#include <tuple> // for tie

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace data {

xt::xtensor<double, 1> compton_profile_pz;

std::vector<PhotonInteraction> elements;
std::unordered_map<std::string, int> element_map;

} // namespace data

//==============================================================================
// PhotonInteraction implementation
//==============================================================================

PhotonInteraction::PhotonInteraction(hid_t group, int i_element)
  : i_element_{i_element}
{
  // Get name of nuclide from group, removing leading '/'
  name_ = object_name(group).substr(1);

  // Get atomic number
  read_attribute(group, "Z", Z_);

  // Determine number of energies and read energy grid
  read_dataset(group, "energy", energy_);

  // Read coherent scattering
  hid_t rgroup = open_group(group, "coherent");
  read_dataset(rgroup, "xs", coherent_);

  hid_t dset = open_dataset(rgroup, "integrated_scattering_factor");
  coherent_int_form_factor_ = Tabulated1D{dset};
  close_dataset(dset);

  if (object_exists(group, "anomalous_real")) {
    dset = open_dataset(rgroup, "anomalous_real");
    coherent_anomalous_real_ = Tabulated1D{dset};
    close_dataset(dset);
  }

  if (object_exists(group, "anomalous_imag")) {
    dset = open_dataset(rgroup, "anomalous_imag");
    coherent_anomalous_imag_ = Tabulated1D{dset};
    close_dataset(dset);
  }
  close_group(rgroup);

  // Read incoherent scattering
  rgroup = open_group(group, "incoherent");
  read_dataset(rgroup, "xs", incoherent_);
  dset = open_dataset(rgroup, "scattering_factor");
  incoherent_form_factor_ = Tabulated1D{dset};
  close_dataset(dset);
  close_group(rgroup);

  // Read pair production
  rgroup = open_group(group, "pair_production_electron");
  read_dataset(rgroup, "xs", pair_production_electron_);
  close_group(rgroup);

  // Read pair production
  if (object_exists(group, "pair_production_nuclear")) {
    rgroup = open_group(group, "pair_production_nuclear");
    read_dataset(rgroup, "xs", pair_production_nuclear_);
    close_group(rgroup);
  } else {
    pair_production_nuclear_ = xt::zeros_like(energy_);
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
    heating_ = xt::zeros_like(energy_);
  }

  // Read subshell photoionization cross section and atomic relaxation data
  rgroup = open_group(group, "subshells");
  std::vector<std::string> designators;
  read_attribute(rgroup, "designators", designators);
  auto n_shell = designators.size();
  if (n_shell == 0) {
    throw std::runtime_error{"Photoatomic data for " + name_ +
      " does not have subshell data."};
  }

  for (int i = 0; i < n_shell; ++i) {
    const auto& designator {designators[i]};

    // TODO: Move to ElectronSubshell constructor

    // Add empty shell
    shells_.emplace_back();

    // Create mapping from designator to index
    int j = 1;
    for (const auto& subshell : SUBSHELLS) {
      if (designator == subshell) {
        shell_map_[j] = i;
        shells_[i].index_subshell = j;
        break;
      }
      ++j;
    }

    // Read binding energy and number of electrons
    auto& shell {shells_.back()};
    hid_t tgroup = open_group(rgroup, designator.c_str());
    read_attribute(tgroup, "binding_energy", shell.binding_energy);
    read_attribute(tgroup, "num_electrons", shell.n_electrons);

    // Read subshell cross section
    dset = open_dataset(tgroup, "xs");
    read_attribute(dset, "threshold_idx", j);
    shell.threshold = j;
    close_dataset(dset);
    read_dataset(tgroup, "xs", shell.cross_section);

    auto& xs = shell.cross_section;
    xs = xt::where(xs > 0.0, xt::log(xs), -500.0);

    if (object_exists(tgroup, "transitions")) {
      // Determine dimensions of transitions
      dset = open_dataset(tgroup, "transitions");
      auto dims = object_shape(dset);
      close_dataset(dset);

      int n_transition = dims[0];
      shell.n_transitions = n_transition;
      if (n_transition > 0) {
        xt::xtensor<double, 2> matrix;
        read_dataset(tgroup, "transitions", matrix);

        shell.transition_subshells = xt::view(matrix, xt::all(), xt::range(0, 2));
        shell.transition_energy = xt::view(matrix, xt::all(), 2);
        shell.transition_probability = xt::view(matrix, xt::all(), 3);
        shell.transition_probability /= xt::sum(shell.transition_probability)();
      }
    } else {
      shell.n_transitions = 0;
    }
    close_group(tgroup);
  }
  close_group(rgroup);

  // Determine number of electron shells
  rgroup = open_group(group, "compton_profiles");

  // Read electron shell PDF and binding energies
  read_dataset(rgroup, "num_electrons", electron_pdf_);
  electron_pdf_ /= xt::sum(electron_pdf_);
  read_dataset(rgroup, "binding_energy", binding_energy_);

  // Read Compton profiles
  read_dataset(rgroup, "J", profile_pdf_);

  // Get Compton profile momentum grid
  if (data::compton_profile_pz.size() == 0) {
    read_dataset(rgroup, "pz", data::compton_profile_pz);
  }
  close_group(rgroup);

  // Create Compton profile CDF
  auto n_profile = data::compton_profile_pz.size();
  profile_cdf_ = xt::empty<double>({n_shell, n_profile});
  for (int i = 0; i < n_shell; ++i) {
    double c = 0.0;
    profile_cdf_(i,0) = 0.0;
    for (int j = 0; j < n_profile - 1; ++j) {
      c += 0.5*(data::compton_profile_pz(j+1) - data::compton_profile_pz(j)) *
        (profile_pdf_(i,j) + profile_pdf_(i,j+1));
      profile_cdf_(i,j+1) = c;
    }
  }

  // Calculate total pair production
  pair_production_total_ = pair_production_nuclear_ + pair_production_electron_;

  if (settings::electron_treatment == ElectronTreatment::TTB) {
    // Read bremsstrahlung scaled DCS
    rgroup = open_group(group, "bremsstrahlung");
    read_dataset(rgroup, "dcs", dcs_);
    auto n_e = dcs_.shape()[0];
    auto n_k = dcs_.shape()[1];

    // Get energy grids used for bremsstrahlung DCS and for stopping powers
    xt::xtensor<double, 1> electron_energy;
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
    int photon = static_cast<int>(Particle::Type::photon);
    const auto& E {electron_energy};
    double cutoff = settings::energy_cutoff[photon];
    if (cutoff > E(0)) {
      size_t i_grid = lower_bound_index(E.cbegin(), E.cend(),
        settings::energy_cutoff[photon]);

      // calculate interpolation factor
      double f = (std::log(cutoff) - std::log(E(i_grid))) /
        (std::log(E(i_grid+1)) - std::log(E(i_grid)));

      // Interpolate bremsstrahlung DCS at the cutoff energy and truncate
      xt::xtensor<double, 2> dcs({n_e - i_grid, n_k});
      for (int i = 0; i < n_k; ++i) {
        double y = std::exp(std::log(dcs_(i_grid,i)) +
              f*(std::log(dcs_(i_grid+1,i)) - std::log(dcs_(i_grid,i))));
        auto col_i = xt::view(dcs, xt::all(), i);
        col_i(0) = y;
        for (int j = i_grid + 1; j < n_e; ++j) {
          col_i(j - i_grid) = dcs_(j, i);
        }
      }
      dcs_ = dcs;

      xt::xtensor<double, 1> frst {cutoff};
      electron_energy = xt::concatenate(xt::xtuple(
        frst, xt::view(electron_energy, xt::range(i_grid+1, n_e))));
    }

    // Set incident particle energy grid
    if (data::ttb_e_grid.size() == 0) {
      data::ttb_e_grid = electron_energy;
    }

    // Calculate the radiative stopping power
    stopping_power_radiative_ = xt::empty<double>({data::ttb_e_grid.size()});
    for (int i = 0; i < data::ttb_e_grid.size(); ++i) {
      // Integrate over reduced photon energy
      double c = 0.0;
      for (int j = 0; j < data::ttb_k_grid.size() - 1; ++j) {
        c += 0.5 * (dcs_(i, j+1) + dcs_(i, j)) * (data::ttb_k_grid(j+1) -
          data::ttb_k_grid(j));
      }
      double e = data::ttb_e_grid(i);

      // Square of the ratio of the speed of light to the velocity of the
      // charged particle
      double beta_sq = e * (e + 2.0 * MASS_ELECTRON_EV) / ((e +
        MASS_ELECTRON_EV) * (e + MASS_ELECTRON_EV));

      stopping_power_radiative_(i) = Z_ * Z_ / beta_sq * e * c;
    }
  }

  // Take logarithm of energies and cross sections since they are log-log
  // interpolated
  energy_ = xt::log(energy_);
  coherent_ = xt::where(coherent_ > 0.0, xt::log(coherent_), -500.0);
  incoherent_ = xt::where(incoherent_ > 0.0, xt::log(incoherent_), -500.0);
  photoelectric_total_ = xt::where(photoelectric_total_ > 0.0,
    xt::log(photoelectric_total_), -500.0);
  pair_production_total_ = xt::where(pair_production_total_ > 0.0,
    xt::log(pair_production_total_), -500.0);
  heating_ = xt::where(heating_ > 0.0, xt::log(heating_), -500.0);
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
    double x = MASS_ELECTRON_EV/PLANCK_C*alpha*std::sqrt(0.5*(1.0 - *mu));

    // Calculate S(x, Z) and S(x_max, Z)
    double form_factor_x = incoherent_form_factor_(x);
    if (form_factor_xmax == 0.0) {
      form_factor_xmax = incoherent_form_factor_(MASS_ELECTRON_EV/PLANCK_C*alpha);
    }

    // Perform rejection on form factor
    if (prn(seed) < form_factor_x / form_factor_xmax) {
      if (doppler) {
        double E_out;
        this->compton_doppler(alpha, *mu, &E_out, i_shell, seed);
        *alpha_out = E_out/MASS_ELECTRON_EV;
      } else {
        *i_shell = -1;
      }
      break;
    }
  }
}

void PhotonInteraction::compton_doppler(double alpha, double mu,
  double* E_out, int* i_shell, uint64_t* seed) const
{
  auto n = data::compton_profile_pz.size();

  int shell; // index for shell
  while (true) {
    // Sample electron shell
    double rn = prn(seed);
    double c = 0.0;
    for (shell = 0; shell < electron_pdf_.size(); ++shell) {
      c += electron_pdf_(shell);
      if (rn < c) break;
    }

    // Determine binding energy of shell
    double E_b = binding_energy_(shell);

    // Determine p_z,max
    double E = alpha*MASS_ELECTRON_EV;
    if (E < E_b) {
      *E_out = alpha/(1 + alpha*(1 - mu))*MASS_ELECTRON_EV;
      break;
    }

    double pz_max = -FINE_STRUCTURE*(E_b - (E - E_b)*alpha*(1.0 - mu)) /
      std::sqrt(2.0*E*(E - E_b)*(1.0 - mu) + E_b*E_b);
    if (pz_max < 0.0) {
      *E_out = alpha/(1 + alpha*(1 - mu))*MASS_ELECTRON_EV;
      break;
    }

    // Determine profile cdf value corresponding to p_z,max
    double c_max;
    if (pz_max > data::compton_profile_pz(n - 1)) {
      c_max = profile_cdf_(shell, n - 1);
    } else {
      int i = lower_bound_index(data::compton_profile_pz.cbegin(),
        data::compton_profile_pz.cend(), pz_max);
      double pz_l = data::compton_profile_pz(i);
      double pz_r = data::compton_profile_pz(i + 1);
      double p_l = profile_pdf_(shell, i);
      double p_r = profile_pdf_(shell, i + 1);
      double c_l = profile_cdf_(shell, i);
      if (pz_l == pz_r) {
        c_max = c_l;
      } else if (p_l == p_r) {
        c_max = c_l + (pz_max - pz_l)*p_l;
      } else {
        double m = (p_l - p_r)/(pz_l - pz_r);
        c_max = c_l + (std::pow((m*(pz_max - pz_l) + p_l), 2) - p_l*p_l)/(2.0*m);
      }
    }

    // Sample value on bounded cdf
    c = prn(seed)*c_max;

    // Determine pz corresponding to sampled cdf value
    auto cdf_shell = xt::view(profile_cdf_, shell, xt::all());
    int i = lower_bound_index(cdf_shell.cbegin(), cdf_shell.cend(), c);
    double pz_l = data::compton_profile_pz(i);
    double pz_r = data::compton_profile_pz(i + 1);
    double p_l = profile_pdf_(shell, i);
    double p_r = profile_pdf_(shell, i + 1);
    double c_l = profile_cdf_(shell, i);
    double pz;
    if (pz_l == pz_r) {
      pz = pz_l;
    } else if (p_l == p_r) {
      pz = pz_l + (c - c_l)/p_l;
    } else {
      double m = (p_l - p_r)/(pz_l - pz_r);
      pz = pz_l + (std::sqrt(p_l*p_l + 2.0*m*(c - c_l)) - p_l)/m;
    }

    // Determine outgoing photon energy corresponding to electron momentum
    // (solve Eq. 39 in LA-UR-04-0487 for E')
    double momentum_sq = std::pow((pz/FINE_STRUCTURE), 2);
    double f = 1.0 + alpha*(1.0 - mu);
    double a = momentum_sq - f*f;
    double b = 2.0*E*(f - momentum_sq*mu);
    c = E*E*(momentum_sq - 1.0);

    double quad = b*b - 4.0*a*c;
    if (quad < 0) {
      *E_out = alpha/(1 + alpha*(1 - mu))*MASS_ELECTRON_EV;
      break;
    }
    quad = std::sqrt(quad);
    double E_out1 = -(b + quad)/(2.0*a);
    double E_out2 = -(b - quad)/(2.0*a);

    // Determine solution to quadratic equation that is positive
    if (E_out1 > 0.0) {
      if (E_out2 > 0.0) {
        // If both are positive, pick one at random
        *E_out = prn(seed) < 0.5 ? E_out1 : E_out2;
      } else {
        *E_out = E_out1;
      }
    } else {
      if (E_out2 > 0.0) {
        *E_out = E_out2;
      } else {
        // No positive solution -- resample
        continue;
      }
    }
    if (*E_out < E - E_b) break;
  }

  *i_shell = shell;
}

void PhotonInteraction::calculate_xs(Particle& p) const
{
  // Perform binary search on the element energy grid in order to determine
  // which points to interpolate between
  int n_grid = energy_.size();
  double log_E = std::log(p.E_);
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
  if (energy_(i_grid) == energy_(i_grid+1)) ++i_grid;

  // calculate interpolation factor
  double f = (log_E - energy_(i_grid)) / (energy_(i_grid+1) - energy_(i_grid));

  auto& xs {p.photon_xs_[i_element_]};
  xs.index_grid = i_grid;
  xs.interp_factor = f;

  // Calculate microscopic coherent cross section
  xs.coherent = std::exp(coherent_(i_grid) +
    f*(coherent_(i_grid+1) - coherent_(i_grid)));

  // Calculate microscopic incoherent cross section
  xs.incoherent = std::exp(incoherent_(i_grid) +
    f*(incoherent_(i_grid+1) - incoherent_(i_grid)));

  // Calculate microscopic photoelectric cross section
  xs.photoelectric = 0.0;
  for (const auto& shell : shells_) {
    // Check threshold of reaction
    int i_start = shell.threshold;
    if (i_grid < i_start) continue;

    // Evaluation subshell photoionization cross section
    xs.photoelectric +=
      std::exp(shell.cross_section(i_grid-i_start) +
      f*(shell.cross_section(i_grid+1-i_start) -
      shell.cross_section(i_grid-i_start)));
  }

  // Calculate microscopic pair production cross section
  xs.pair_production = std::exp(
    pair_production_total_(i_grid) + f*(
    pair_production_total_(i_grid+1) -
    pair_production_total_(i_grid)));

  // Calculate microscopic total cross section
  xs.total = xs.coherent + xs.incoherent + xs.photoelectric + xs.pair_production;
  xs.last_E = p.E_;
}

double PhotonInteraction::rayleigh_scatter(double alpha, uint64_t* seed) const
{
  double mu;
  while (true) {
    // Determine maximum value of x^2
    double x2_max = std::pow(MASS_ELECTRON_EV/PLANCK_C*alpha, 2);

    // Determine F(x^2_max, Z)
    double F_max = coherent_int_form_factor_(x2_max);

    // Sample cumulative distribution
    double F = prn(seed)*F_max;

    // Determine x^2 corresponding to F
    const auto& x {coherent_int_form_factor_.x()};
    const auto& y {coherent_int_form_factor_.y()};
    int i = lower_bound_index(y.cbegin(), y.cend(), F);
    double r = (F - y[i]) / (y[i+1] - y[i]);
    double x2 = x[i] + r*(x[i+1] - x[i]);

    // Calculate mu
    mu = 1.0 - 2.0*x2/x2_max;

    if (prn(seed) < 0.5*(1.0 + mu*mu)) break;
  }
  return mu;
}

void PhotonInteraction::pair_production(double alpha, double* E_electron,
  double* E_positron, double* mu_electron, double* mu_positron,
  uint64_t* seed) const
{
  constexpr double r[] {
    122.81, 73.167, 69.228, 67.301, 64.696, 61.228,
    57.524, 54.033, 50.787, 47.851, 46.373, 45.401,
    44.503, 43.815, 43.074, 42.321, 41.586, 40.953,
    40.524, 40.256, 39.756, 39.144, 38.462, 37.778,
    37.174, 36.663, 35.986, 35.317, 34.688, 34.197,
    33.786, 33.422, 33.068, 32.740, 32.438, 32.143,
    31.884, 31.622, 31.438, 31.142, 30.950, 30.758,
    30.561, 30.285, 30.097, 29.832, 29.581, 29.411,
    29.247, 29.085, 28.930, 28.721, 28.580, 28.442,
    28.312, 28.139, 27.973, 27.819, 27.675, 27.496,
    27.285, 27.093, 26.911, 26.705, 26.516, 26.304,
    26.108, 25.929, 25.730, 25.577, 25.403, 25.245,
    25.100, 24.941, 24.790, 24.655, 24.506, 24.391,
    24.262, 24.145, 24.039, 23.922, 23.813, 23.712,
    23.621, 23.523, 23.430, 23.331, 23.238, 23.139,
    23.048, 22.967, 22.833, 22.694, 22.624, 22.545,
    22.446, 22.358, 22.264};

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
  double c = a*a*(1.0/(1.0 + a*a) + 0.202059 + a*a*(-0.03693 + a*a*(0.00835 +
    a*a*(-0.00201 + a*a*(0.00049 +  a*a*(-0.00012 + a*a*0.00003))))));

  // The analytical approximation of the DCS underestimates the cross section
  // at low energies. The correction factor f compensates for this.
  double q = std::sqrt(2.0/alpha);
  double f = q*(-0.1774 - 12.10*a + 11.18*a*a)
    + q*q*(8.523 + 73.26*a - 44.41*a*a)
    + q*q*q*(-13.52 - 121.1*a + 96.41*a*a)
    + q*q*q*q*(8.946 + 62.05*a - 63.41*a*a);

  // Calculate phi_1(1/2) and phi_2(1/2). The unnormalized PDF for the reduced
  // energy is given by p = 2*(1/2 - e)^2*phi_1(e) + phi_2(e), where phi_1 and
  // phi_2 are non-negative and maximum at e = 1/2.
  double b = 2.0*r[Z_]/alpha;
  double t1 = 2.0*std::log(1.0 + b*b);
  double t2 = b*std::atan(1.0/b);
  double t3 = b*b*(4.0 - 4.0*t2 - 3.0*std::log(1.0 + 1.0/(b*b)));
  double t4 = 4.0*std::log(r[Z_]) - 4.0*c + f;
  double phi1_max = 7.0/3.0 - t1 - 6.0*t2 - t3 + t4;
  double phi2_max = 11.0/6.0 - t1 - 3.0*t2 + 0.5*t3 + t4;

  // To aid sampling, the unnormalized PDF can be expressed as
  // p = u_1*U_1(e)*pi_1(e) + u_2*U_2(e)*pi_2(e), where pi_1 and pi_2 are
  // normalized PDFs on the interval (e_min, e_max) from which values of e can
  // be sampled using the inverse transform method, and
  // U_1 = phi_1(e)/phi_1(1/2) and U_2 = phi_2(e)/phi_2(1/2) are valid
  // rejection functions. The reduced energy can now be sampled using a
  // combination of the composition and rejection methods.
  double u1 = 2.0/3.0*std::pow(0.5 - 1.0/alpha, 2)*phi1_max;
  double u2 = phi2_max;
  double e;
  while (true) {
    double rn = prn(seed);

    // Sample the index i in (1, 2) using the point probabilities
    // p(1) = u_1/(u_1 + u_2) and p(2) = u_2/(u_1 + u_2)
    int i;
    if (prn(seed) < u1/(u1 + u2)) {
      i = 1;

      // Sample e from pi_1 using the inverse transform method
      e = rn >= 0.5 ?
        0.5 + (0.5 - 1.0/alpha)*std::pow(2.0*rn - 1.0, 1.0/3.0) :
        0.5 - (0.5 - 1.0/alpha)*std::pow(1.0 - 2.0*rn, 1.0/3.0);
    } else {
      i = 2;

      // Sample e from pi_2 using the inverse transform method
      e = 1.0/alpha + (0.5 - 1.0/alpha)*2.0*rn;
    }

    // Calculate phi_i(e) and deliver e if rn <= U_i(e)
    b = r[Z_]/(2.0*alpha*e*(1.0 - e));
    t1 = 2.0*std::log(1.0 + b*b);
    t2 = b*std::atan(1.0/b);
    t3 = b*b*(4.0 - 4.0*t2 - 3.0*std::log(1.0 + 1.0/(b*b)));
    if (i == 1) {
      double phi1 = 7.0/3.0 - t1 - 6.0*t2 - t3 + t4;
      if (prn(seed) <= phi1/phi1_max) break;
    } else {
      double phi2 = 11.0/6.0 - t1 - 3.0*t2 + 0.5*t3 + t4;
      if (prn(seed) <= phi2/phi2_max) break;
    }
  }

  // Compute the kinetic energy of the electron and the positron
  *E_electron = (alpha*e - 1.0)*MASS_ELECTRON_EV;
  *E_positron = (alpha*(1.0 - e) - 1.0)*MASS_ELECTRON_EV;

  // Sample the scattering angle of the electron. The cosine of the polar
  // angle of the direction relative to the incident photon is sampled from
  // p(mu) = C/(1 - beta*mu)^2 using the inverse transform method.
  double beta = std::sqrt(*E_electron*(*E_electron + 2.0*MASS_ELECTRON_EV))
    / (*E_electron + MASS_ELECTRON_EV)  ;
  double rn = 2.0*prn(seed) - 1.0;
  *mu_electron = (rn + beta)/(rn*beta + 1.0);

  // Sample the scattering angle of the positron
  beta = std::sqrt(*E_positron*(*E_positron + 2.0*MASS_ELECTRON_EV))
    / (*E_positron + MASS_ELECTRON_EV);
  rn = 2.0*prn(seed) - 1.0;
  *mu_positron = (rn + beta)/(rn*beta + 1.0);
}

void PhotonInteraction::atomic_relaxation(const ElectronSubshell& shell, Particle& p) const
{
  // If no transitions, assume fluorescent photon from captured free electron
  if (shell.n_transitions == 0) {
    double mu = 2.0*prn(p.current_seed()) - 1.0;
    double phi = 2.0*PI*prn(p.current_seed());
    Direction u;
    u.x = mu;
    u.y = std::sqrt(1.0 - mu*mu)*std::cos(phi);
    u.z = std::sqrt(1.0 - mu*mu)*std::sin(phi);
    double E = shell.binding_energy;
    p.create_secondary(u, E, Particle::Type::photon);
    return;
  }

  // Sample transition
  double rn = prn(p.current_seed());
  double c = 0.0;
  int i_transition;
  for (i_transition = 0; i_transition < shell.n_transitions; ++i_transition) {
    c += shell.transition_probability(i_transition);
    if (rn < c) break;
  }

  // Get primary and secondary subshell designators
  int primary = shell.transition_subshells(i_transition, 0);
  int secondary = shell.transition_subshells(i_transition, 1);

  // Sample angle isotropically
  double mu = 2.0*prn(p.current_seed()) - 1.0;
  double phi = 2.0*PI*prn(p.current_seed());
  Direction u;
  u.x = mu;
  u.y = std::sqrt(1.0 - mu*mu)*std::cos(phi);
  u.z = std::sqrt(1.0 - mu*mu)*std::sin(phi);

  // Get the transition energy
  double E = shell.transition_energy(i_transition);

  if (secondary != 0) {
    // Non-radiative transition -- Auger/Coster-Kronig effect

    // Create auger electron
    p.create_secondary(u, E, Particle::Type::electron);

    // Fill hole left by emitted auger electron
    int i_hole = shell_map_.at(secondary);
    const auto& hole = shells_[i_hole];
    this->atomic_relaxation(hole, p);
  } else {
    // Radiative transition -- get X-ray energy

    // Create fluorescent photon
    p.create_secondary(u, E, Particle::Type::photon);
  }

  // Fill hole created by electron transitioning to the photoelectron hole
  int i_hole = shell_map_.at(primary);
  const auto& hole = shells_[i_hole];
  this->atomic_relaxation(hole, p);
}

//==============================================================================
// Non-member functions
//==============================================================================

std::pair<double, double> klein_nishina(double alpha, uint64_t* seed)
{
  double alpha_out, mu;
  double beta = 1.0 + 2.0*alpha;
  if (alpha < 3.0) {
    // Kahn's rejection method
    double t = beta/(beta + 8.0);
    double x;
    while (true) {
      if (prn(seed) < t) {
        // Left branch of flow chart
        double r = 2.0*prn(seed);
        x = 1.0 + alpha*r;
        if (prn(seed) < 4.0/x*(1.0 - 1.0/x)) {
          mu = 1 - r;
          break;
        }
      } else {
        // Right branch of flow chart
        x = beta/(1.0 + 2.0*alpha*prn(seed));
        mu = 1.0 + (1.0 - x)/alpha;
        if (prn(seed) < 0.5*(mu*mu + 1.0/x)) break;
      }
    }
    alpha_out = alpha/x;

  } else {
    // Koblinger's direct method
    double gamma = 1.0 - std::pow(beta, -2);
    double s = prn(seed)*(4.0/alpha + 0.5*gamma +
      (1.0 - (1.0 + beta)/(alpha*alpha))*std::log(beta));
    if (s <= 2.0/alpha) {
      // For first term, x = 1 + 2ar
      // Therefore, a' = a/(1 + 2ar)
      alpha_out = alpha/(1.0 + 2.0*alpha*prn(seed));
    } else if (s <= 4.0/alpha) {
      // For third term, x = beta/(1 + 2ar)
      // Therefore, a' = a(1 + 2ar)/beta
      alpha_out = alpha*(1.0 + 2.0*alpha*prn(seed))/beta;
    } else if (s <= 4.0/alpha + 0.5*gamma) {
      // For fourth term, x = 1/sqrt(1 - gamma*r)
      // Therefore, a' = a*sqrt(1 - gamma*r)
      alpha_out = alpha*std::sqrt(1.0 - gamma*prn(seed));
    } else {
      // For third term, x = beta^r
      // Therefore, a' = a/beta^r
      alpha_out = alpha/std::pow(beta, prn(seed));
    }

    // Calculate cosine of scattering angle based on basic relation
    mu = 1.0 + 1.0/alpha - 1.0/alpha_out;
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
