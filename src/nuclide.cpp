#include "openmc/nuclide.h"

#include "openmc/container_util.h"
#include "openmc/endf.h"
#include "openmc/error.h"
#include "openmc/hdf5_interface.h"
#include "openmc/message_passing.h"
#include "openmc/random_lcg.h"
#include "openmc/search.h"
#include "openmc/settings.h"
#include "openmc/string_utils.h"

#include <algorithm> // for sort
#include <string> // for to_string, stoi

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace data {
std::array<double, 2> energy_min {0.0, 0.0};
std::array<double, 2> energy_max {INFTY, INFTY};
std::vector<std::unique_ptr<Nuclide>> nuclides;
} // namespace data

namespace simulation {
NuclideMicroXS* micro_xs;
MaterialMacroXS material_xs;
} // namespace simulation

//==============================================================================
// Nuclide implementation
//==============================================================================

Nuclide::Nuclide(hid_t group, const double* temperature, int n, int i_nuclide)
{
  // Get name of nuclide from group, removing leading '/'
  name_ = object_name(group).substr(1);

  read_attribute(group, "Z", Z_);
  read_attribute(group, "A", A_);
  read_attribute(group, "metastable", metastable_);
  read_attribute(group, "atomic_weight_ratio", awr_);
  i_nuclide_ = i_nuclide;

  // Determine temperatures available
  hid_t kT_group = open_group(group, "kTs");
  auto dset_names = dataset_names(kT_group);
  std::vector<double> temps_available;
  for (const auto& name : dset_names) {
    double T;
    read_dataset(kT_group, name.c_str(), T);
    temps_available.push_back(T / K_BOLTZMANN);
  }
  std::sort(temps_available.begin(), temps_available.end());

  // If only one temperature is available, revert to nearest temperature
  if (temps_available.size() == 1 && settings::temperature_method == TEMPERATURE_INTERPOLATION) {
    if (mpi::master) {
      warning("Cross sections for " + name_ + " are only available at one "
        "temperature. Reverting to nearest temperature method.");
    }
    settings::temperature_method = TEMPERATURE_NEAREST;
  }

  // Determine actual temperatures to read -- start by checking whether a
  // temperature range was given, in which case all temperatures in the range
  // are loaded irrespective of what temperatures actually appear in the model
  std::vector<int> temps_to_read;
  double T_min = n > 0 ? settings::temperature_range[0] : 0.0;
  double T_max = n > 0 ? settings::temperature_range[1] : INFTY;
  if (T_max > 0.0) {
    for (auto T : temps_available) {
      if (T_min <= T && T <= T_max) {
        temps_to_read.push_back(std::round(T));
      }
    }
  }

  switch (settings::temperature_method) {
  case TEMPERATURE_NEAREST:
    // Find nearest temperatures
    for (int i = 0; i < n; ++i) {
      double T_desired = temperature[i];

      // Determine closest temperature
      double min_delta_T = INFTY;
      double T_actual;
      for (auto T : temps_available) {
        double delta_T = std::abs(T - T_desired);
        if (delta_T < min_delta_T) {
          T_actual = T;
          min_delta_T = delta_T;
        }
      }

      if (std::abs(T_actual - T_desired) < settings::temperature_tolerance) {
        if (!contains(temps_to_read, std::round(T_actual))) {
          temps_to_read.push_back(std::round(T_actual));

          // Write warning for resonance scattering data if 0K is not available
          if (std::abs(T_actual - T_desired) > 0 && T_desired == 0 && mpi::master) {
            warning(name_ + " does not contain 0K data needed for resonance "
              "scattering options selected. Using data at " + std::to_string(T_actual)
              + " K instead.");
          }
        }
      } else {
        fatal_error("Nuclear data library does not contain cross sections for " +
          name_ + " at or near " + std::to_string(T_desired) + " K.");
      }
    }
    break;

  case TEMPERATURE_INTERPOLATION:
    // If temperature interpolation or multipole is selected, get a list of
    // bounding temperatures for each actual temperature present in the model
    for (int i = 0; i < n; ++i) {
      double T_desired = temperature[i];

      bool found_pair = false;
      for (int j = 0; j < temps_available.size() - 1; ++j) {
        if (temps_available[j] <= T_desired && T_desired < temps_available[j + 1]) {
          int T_j = std::round(temps_available[j]);
          int T_j1 = std::round(temps_available[j+1]);
          if (!contains(temps_to_read, T_j)) {
            temps_to_read.push_back(T_j);
          }
          if (!contains(temps_to_read, T_j1)) {
            temps_to_read.push_back(T_j1);
          }
          found_pair = true;
        }
      }

      if (!found_pair) {
        fatal_error("Nuclear data library does not contain cross sections for " +
          name_ +" at temperatures that bound " + std::to_string(T_desired) + " K.");
      }
    }
    break;
  }

  // Sort temperatures to read
  std::sort(temps_to_read.begin(), temps_to_read.end());

  hid_t energy_group = open_group(group, "energy");
  for (const auto& T : temps_to_read) {
    std::string dset {std::to_string(T) + "K"};

    // Determine exact kT values
    double kT;
    read_dataset(kT_group, dset.c_str(), kT);
    kTs_.push_back(kT);

    // Read energy grid
    grid_.emplace_back();
    read_dataset(energy_group, dset.c_str(), grid_.back().energy);
  }
  close_group(kT_group);

  // Check for 0K energy grid
  if (object_exists(energy_group, "0K")) {
    read_dataset(energy_group, "0K", energy_0K_);
  }
  close_group(energy_group);

  // Read reactions
  hid_t rxs_group = open_group(group, "reactions");
  for (auto name : group_names(rxs_group)) {
    if (starts_with(name, "reaction_")) {
      hid_t rx_group = open_group(rxs_group, name.c_str());
      reactions_.push_back(std::make_unique<Reaction>(rx_group, temps_to_read));

      // Check for 0K elastic scattering
      const auto& rx = reactions_.back();
      if (rx->mt_ == ELASTIC) {
        if (object_exists(rx_group, "0K")) {
          hid_t temp_group = open_group(rx_group, "0K");
          read_dataset(temp_group, "xs", elastic_0K_);
          close_group(temp_group);
        }
      }
      close_group(rx_group);

      // Determine reaction indices for inelastic scattering reactions
      if (is_inelastic_scatter(rx->mt_) && !rx->redundant_) {
        index_inelastic_scatter_.push_back(reactions_.size() - 1);
      }
    }
  }
  close_group(rxs_group);

  // Read unresolved resonance probability tables if present
  if (object_exists(group, "urr")) {
    urr_present_ = true;
    urr_data_.resize(temps_to_read.size());

    for (int i = 0; i < temps_to_read.size(); i++) {
      // Get temperature as a string
      std::string temp_str {std::to_string(temps_to_read[i]) + "K"};

      // Read probability tables for i-th temperature
      hid_t urr_group = open_group(group, ("urr/" + temp_str).c_str());
      urr_data_[i].from_hdf5(urr_group);
      close_group(urr_group);

      // Check for negative values
      if (xt::any(urr_data_[i].prob_ < 0.) && mpi::master) {
        warning("Negative value(s) found on probability table for nuclide " +
                name_ + " at " + temp_str);
      }
    }

    // If the inelastic competition flag indicates that the inelastic cross
    // section should be determined from a normal reaction cross section, we
    // need to get the index of the reaction.
    if (temps_to_read.size() > 0) {
      if (urr_data_[0].inelastic_flag_ > 0) {
        for (int i = 0; i < reactions_.size(); i++) {
          if (reactions_[i]->mt_ == urr_data_[0].inelastic_flag_) {
            urr_inelastic_ = i;
          }
        }

        // Abort if no corresponding inelastic reaction was found
        if (urr_inelastic_ == C_NONE) {
          fatal_error("Could no find inelastic reaction specified on "
                      "unresolved resonance probability table.");
        }
      }
    }
  }

  // Check for nu-total
  if (object_exists(group, "total_nu")) {
    // Read total nu data
    hid_t nu_group = open_group(group, "total_nu");
    hid_t nu_dset = open_dataset(nu_group, "yield");
    std::string func_type;
    read_attribute(nu_dset, "type", func_type);
    if (func_type == "Tabulated1D") {
      total_nu_ = std::make_unique<Tabulated1D>(nu_dset);
    } else if (func_type == "Polynomial") {
      total_nu_ = std::make_unique<Polynomial>(nu_dset);
    }
    close_dataset(nu_dset);
    close_group(nu_group);
  }

  this->create_derived();
}

void Nuclide::create_derived()
{
  for (int i = 0; i < reactions_.size(); ++i) {
    const auto& rx {reactions_[i]};

    for (int t = 0; t < kTs_.size(); ++t) {
      // Skip redundant reactions
      if (rx->redundant_) continue;

      if (is_fission(rx->mt_)) {
        fissionable_ = true;

        // Keep track of fission reactions
        if (t == 0) {
          fission_rx_.push_back(rx.get());
          if (rx->mt_ == N_F) has_partial_fission_ = true;
        }
      }
    }
  }

  // Determine number of delayed neutron precursors
  if (fissionable_) {
    for (const auto& product : fission_rx_[0]->products_) {
      if (product.emission_mode_ == EmissionMode::delayed) {
        ++n_precursor_;
      }
    }
  }

  if (settings::res_scat_on) {
    // Determine if this nuclide should be treated as a resonant scatterer
    if (!settings::res_scat_nuclides.empty()) {
      // If resonant nuclides were specified, check the list explicitly
      for (const auto& name : settings::res_scat_nuclides) {
        if (name_ == name) {
          resonant_ = true;

          // Make sure nuclide has 0K data
          if (energy_0K_.empty()) {
            fatal_error("Cannot treat " + name_ + " as a resonant scatterer "
              "because 0 K elastic scattering data is not present.");
          }
          break;
        }
      }
    } else {
      // Otherwise, assume that any that have 0 K elastic scattering data are
      // resonant
      resonant_ = !energy_0K_.empty();
    }

    if (resonant_) {
      // Build CDF for 0K elastic scattering
      double xs_cdf_sum = 0.0;
      xs_cdf_.resize(energy_0K_.size());
      xs_cdf_[0] = 0.0;

      const auto& E = energy_0K_;
      auto& xs = elastic_0K_;
      for (int i = 0; i < E.size() - 1; ++i) {
        // Negative cross sections result in a CDF that is not monotonically
        // increasing. Set all negative xs values to zero.
        if (xs[i] < 0.0) xs[i] = 0.0;

        // build xs cdf
        xs_cdf_sum += (std::sqrt(E[i])*xs[i] + std::sqrt(E[i+1])*xs[i+1])
              / 2.0 * (E[i+1] - E[i]);
        xs_cdf_[i] = xs_cdf_sum;
      }
    }
  }
}

double Nuclide::nu(double E, EmissionMode mode, int group) const
{
  if (!fissionable_) return 0.0;

  switch (mode) {
  case EmissionMode::prompt:
    return (*fission_rx_[0]->products_[0].yield_)(E);
  case EmissionMode::delayed:
    if (n_precursor_ > 0) {
      auto rx = fission_rx_[0];
      if (group >= 1 && group < rx->products_.size()) {
        // If delayed group specified, determine yield immediately
        return (*rx->products_[group].yield_)(E);
      } else {
        double nu {0.0};

        for (int i = 1; i < rx->products_.size(); ++i) {
          // Skip any non-neutron products
          const auto& product = rx->products_[i];
          if (product.particle_ != ParticleType::neutron) continue;

          // Evaluate yield
          if (product.emission_mode_ == EmissionMode::delayed) {
            nu += (*product.yield_)(E);
          }
        }
        return nu;
      }
    } else {
      return 0.0;
    }
  case EmissionMode::total:
    if (total_nu_) {
      return (*total_nu_)(E);
    } else {
      return (*fission_rx_[0]->products_[0].yield_)(E);
    }
  }
}

void Nuclide::calculate_elastic_xs() const
{
  // Get temperature index, grid index, and interpolation factor
  auto& micro = simulation::micro_xs[i_nuclide_];
  int i_temp = micro.index_temp - 1;
  int i_grid = micro.index_grid - 1;
  double f = micro.interp_factor;

  if (i_temp >= 0) {
    const auto& xs = reactions_[0]->xs_[i_temp].value;
    micro.elastic = (1.0 - f)*xs[i_grid] + f*xs[i_grid + 1];
  }
}

double Nuclide::elastic_xs_0K(double E) const
{
  // Determine index on nuclide energy grid
  int i_grid;
  if (E < energy_0K_.front()) {
    i_grid = 0;
  } else if (E > energy_0K_.back()) {
    i_grid = energy_0K_.size() - 2;
  } else {
    i_grid = lower_bound_index(energy_0K_.begin(), energy_0K_.end(), E);
  }

  // check for rare case where two energy points are the same
  if (energy_0K_[i_grid] == energy_0K_[i_grid+1]) ++i_grid;

  // calculate interpolation factor
  double f = (E - energy_0K_[i_grid]) /
    (energy_0K_[i_grid + 1] - energy_0K_[i_grid]);

  // Calculate microscopic nuclide elastic cross section
  return (1.0 - f)*elastic_0K_[i_grid] + f*elastic_0K_[i_grid + 1];
}

void Nuclide::calculate_urr_xs(const int i_temp, const double E)
{
  auto& micro = simulation::micro_xs[i_nuclide_];
  micro.use_ptable = true;

  // Create a shorthand for the URR data
  UrrData* urr = &(urr_data_[i_temp]);

  // Determine the energy table
  int i_energy = 0;
  while (true) {
    if (E < urr->energy_(i_energy + 1)) {break;}
    i_energy++;
  }

  // Sample the probability table using the cumulative distribution

  // Random nmbers for the xs calculation are sampled from a separate stream.
  // This guarantees the rnadomness and, at the same time, makes sure we
  // reuse random numbers for the same nuclide at different temperatures,
  // therefore preserving correlation of temperature in probability tables.
  prn_set_stream(STREAM_URR_PTABLE);
  //TODO: to maintain the same random number stream as the Fortran code this
  //replaces, the seed is set with i_nuclide_ + 1 instead of i_nuclide_
  double r = future_prn(static_cast<int64_t>(i_nuclide_ + 1));
  prn_set_stream(STREAM_TRACKING);

  int i_low = 0;
  while (true) {
    if (urr->prob_(i_energy, URR_CUM_PROB, i_low) > r) {break;}
    i_low++;
  }
  int i_up = 0;
  while (true) {
    if (urr->prob_(i_energy + 1, URR_CUM_PROB, i_up) > r) {break;}
    i_up++;
  }

  // Determine elastic, fission, and capture cross sections from the
  // probability table
  double elastic = 0.;
  double fission = 0.;
  double capture = 0.;
  double f;
  if (urr->interp_ == Interpolation::lin_lin) {
    // Determine the interpolation factor on the table
    f = (E - urr->energy_(i_energy)) /
         (urr->energy_(i_energy + 1) - urr->energy_(i_energy));

    elastic = (1. - f) * urr->prob_(i_energy, URR_ELASTIC, i_low) +
         f * urr->prob_(i_energy + 1, URR_ELASTIC, i_up);
    fission = (1. - f) * urr->prob_(i_energy, URR_FISSION, i_low) +
         f * urr->prob_(i_energy + 1, URR_FISSION, i_up);
    capture = (1. - f) * urr->prob_(i_energy, URR_N_GAMMA, i_low) +
         f * urr->prob_(i_energy + 1, URR_N_GAMMA, i_up);
  } else if (urr->interp_ == Interpolation::log_log) {
    // Determine interpolation factor on the table
    f = std::log(E / urr->energy_(i_energy)) /
         std::log(urr->energy_(i_energy + 1) / urr->energy_(i_energy));

    // Calculate the elastic cross section/factor
    if ((urr->prob_(i_energy, URR_ELASTIC, i_low) > 0.) &&
        (urr->prob_(i_energy + 1, URR_ELASTIC, i_up) > 0.)) {
      elastic =
           std::exp((1. - f) *
                    std::log(urr->prob_(i_energy, URR_ELASTIC, i_low)) +
                    f * std::log(urr->prob_(i_energy + 1, URR_ELASTIC, i_up)));
    } else {
      elastic = 0.;
    }

    // Calculate the fission cross section/factor
    if ((urr->prob_(i_energy, URR_FISSION, i_low) > 0.) &&
        (urr->prob_(i_energy + 1, URR_FISSION, i_up) > 0.)) {
      fission =
           std::exp((1. - f) *
                    std::log(urr->prob_(i_energy, URR_FISSION, i_low)) +
                    f * std::log(urr->prob_(i_energy + 1, URR_FISSION, i_up)));
    } else {
      fission = 0.;
    }

    // Calculate the capture cross section/factor
    if ((urr->prob_(i_energy, URR_N_GAMMA, i_low) > 0.) &&
        (urr->prob_(i_energy + 1, URR_N_GAMMA, i_up) > 0.)) {
      capture =
           std::exp((1. - f) *
                    std::log(urr->prob_(i_energy, URR_N_GAMMA, i_low)) +
                    f * std::log(urr->prob_(i_energy + 1, URR_N_GAMMA, i_up)));
    } else {
      capture = 0.;
    }
  }

  // Determine the treatment of inelastic scattering
  double inelastic = 0.;
  if (urr->inelastic_flag_ != C_NONE) {
    // get interpolation factor
    f = micro.interp_factor;

    // Determine inelastic scattering cross section
    Reaction* rx = reactions_[urr_inelastic_].get();
    int xs_index = micro.index_grid - rx->xs_[i_temp].threshold;
    if (xs_index >= 0) {
      inelastic = (1. - f) * rx->xs_[i_temp].value[xs_index] +
           f * rx->xs_[i_temp].value[xs_index + 1];
    }
  }

  // Multiply by smooth cross-section if needed
  if (urr->multiply_smooth_) {
    calculate_elastic_xs();
    elastic *= micro.elastic;
    capture *= (micro.absorption - micro.fission);
    fission *= micro.fission;
  }

  // Check for negative values
  if (elastic < 0.) {elastic = 0.;}
  if (fission < 0.) {fission = 0.;}
  if (capture < 0.) {capture = 0.;}

  // Set elastic, absorption, fission, and total x/s. Note that the total x/s
  // is calculated as a sum of partials instead of the table-provided value
  micro.elastic = elastic;
  micro.absorption = capture + fission;
  micro.fission = fission;
  micro.total = elastic + inelastic + capture + fission;

  // Determine nu-fission cross-section
  if (fissionable_) {
    micro.nu_fission = nu(E, EmissionMode::total) * micro.fission;
  }

}

//==============================================================================
// Fortran compatibility functions
//==============================================================================

extern "C" void
set_particle_energy_bounds(int particle, double E_min, double E_max)
{
  data::energy_min[particle - 1] = E_min;
  data::energy_max[particle - 1] = E_max;
}

extern "C" Nuclide* nuclide_from_hdf5_c(hid_t group, const double* temperature, int n)
{
  data::nuclides.push_back(std::make_unique<Nuclide>(group, temperature, n,
                                                     data::nuclides.size()));
  return data::nuclides.back().get();
}

extern "C" Reaction* nuclide_reaction(Nuclide* nuc, int i_rx)
{
  return nuc->reactions_[i_rx-1].get();
}

extern "C" void nuclides_clear() { data::nuclides.clear(); }


extern "C" NuclideMicroXS* micro_xs_ptr();
void set_micro_xs()
{
#pragma omp parallel
  {
    simulation::micro_xs = micro_xs_ptr();
  }
}

extern "C" void
nuclide_calculate_urr_xs(const int i_nuclide, const int i_temp, const double E)
{
  data::nuclides[i_nuclide - 1]->calculate_urr_xs(i_temp - 1, E);
}

} // namespace openmc
