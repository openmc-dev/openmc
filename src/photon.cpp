#include "openmc/photon.h"

#include "openmc/constants.h"
#include "openmc/hdf5_interface.h"
#include "openmc/particle.h"
#include "openmc/search.h"
#include "openmc/settings.h"

#include "xtensor/xbuilder.hpp"
#include "xtensor/xoperation.hpp"
#include "xtensor/xview.hpp"

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace data {

xt::xtensor<double, 1> compton_profile_pz;
xt::xtensor<double, 1> ttb_e_grid;
xt::xtensor<double, 1> ttb_k_grid;

std::vector<PhotonInteraction> elements;

} // namespace data

namespace simulation {
ElementMicroXS* micro_photon_xs;
} // namespace simulation

//==============================================================================
// PhotonInteraction implementation
//==============================================================================

PhotonInteraction::PhotonInteraction(hid_t group)
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

  // Read subshell photoionization cross section and atomic relaxation data
  rgroup = open_group(group, "subshells");
  std::vector<std::string> designators;
  read_attribute(rgroup, "designators", designators);
  auto n_shell = designators.size();
  for (int i = 0; i < n_shell; ++i) {
    const auto& designator {designators[i]};

    // TODO: Move to ElectronSubshell constructor

    // Add empty shell
    shells_.emplace_back();

    // Create mapping from designator to index
    int j = 0;
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
    // TODO: off-by-one
    shell.threshold = j - 1;
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

  // Get Compton profile momentum grid. By deafult, an xtensor has a size of 1.
  // TODO: Update version of xtensor and change to 0
  if (data::compton_profile_pz.size() == 1) {
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

  if (settings::electron_treatment == ELECTRON_TTB) {
    // Read bremsstrahlung scaled DCS
    rgroup = open_group(group, "bremsstrahlung");
    read_dataset(rgroup, "dcs", dcs_);
    auto n_e = dcs_.shape()[0];
    auto n_k = dcs_.shape()[1];

    // Get energy grids used for bremsstrahlung DCS and for stopping powers
    xt::xtensor<double, 1> electron_energy;
    read_dataset(rgroup, "electron_energy", electron_energy);
    if (data::ttb_k_grid.size() == 1) {
      read_dataset(rgroup, "photon_energy", data::ttb_k_grid);
    }
    close_group(rgroup);

    // Read stopping power data
    if (Z_ < 99) {
      rgroup = open_group(group, "stopping_powers");
      read_dataset(rgroup, "s_collision", stopping_power_collision_);
      read_dataset(rgroup, "s_radiative", stopping_power_radiative_);
      read_attribute(rgroup, "I", I_);
      close_group(rgroup);
    }

    // Truncate the bremsstrahlung data at the cutoff energy
    int photon = static_cast<int>(ParticleType::photon) - 1;
    const auto& E {electron_energy};
    double cutoff = settings::energy_cutoff[photon];
    if (cutoff > E(0)) {
      size_t i_grid = lower_bound_index(E.cbegin(), E.cend(),
        settings::energy_cutoff[photon]);

      // calculate interpolation factor
      double f = (std::log(cutoff) - std::log(E(i_grid))) /
        (std::log(E(i_grid+1)) - std::log(E(i_grid)));

      // Interpolate collision stopping power at the cutoff energy and truncate
      auto& s_col {stopping_power_collision_};
      double y = std::exp(std::log(s_col(i_grid)) + f*(std::log(s_col(i_grid+1)) -
        std::log(s_col(i_grid))));
      xt::xtensor<double, 1> frst {y};
      stopping_power_collision_ = xt::concatenate(xt::xtuple(
        frst, xt::view(s_col, xt::range(i_grid+1, n_e))));

      // Interpolate radiative stopping power at the cutoff energy and truncate
      auto& s_rad {stopping_power_radiative_};
      y = std::exp(std::log(s_rad(i_grid)) + f*(std::log(s_rad(i_grid+1)) -
        std::log(s_rad(i_grid))));
      frst(0) = y;
      stopping_power_radiative_ = xt::concatenate(xt::xtuple(
        frst, xt::view(s_rad, xt::range(i_grid+1, n_e))));

      // Interpolate bremsstrahlung DCS at the cutoff energy and truncate
      xt::xtensor<double, 2> dcs = xt::empty<double>({n_e - i_grid, n_k});
      for (int i = 0; i < n_k; ++i) {
        y = std::exp(std::log(dcs_(i_grid,i)) +
              f*(std::log(dcs_(i_grid+1,i)) - std::log(dcs_(i_grid,i))));
        auto col_i = xt::view(dcs, xt::all(), i);
        col_i(0) = y;
        for (int j = i_grid + 1; j < n_e; ++j) {
          col_i(j - i_grid) = dcs_(j, i);
        }
      }
      dcs_ = dcs;

      frst(0) = cutoff;
      electron_energy = xt::concatenate(xt::xtuple(
        frst, xt::view(electron_energy, xt::range(i_grid+1, n_e))));
    }

    // Set incident particle energy grid
    if (data::ttb_e_grid.size() == 1) {
      data::ttb_e_grid = electron_energy;
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
}

//==============================================================================
// Fortran compatibility
//==============================================================================

extern "C" void photon_from_hdf5_c(hid_t group)
{
  data::elements.emplace_back(group);
}

} // namespace openmc
