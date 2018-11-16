#include "openmc/nuclide.h"

#include "openmc/container_util.h"
#include "openmc/error.h"
#include "openmc/hdf5_interface.h"
#include "openmc/message_passing.h"
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

//==============================================================================
// Nuclide implementation
//==============================================================================

Nuclide::Nuclide(hid_t group, const double* temperature, int n)
{
  // Get name of nuclide from group, removing leading '/'
  name_ = object_name(group).substr(1);

  read_attribute(group, "Z", Z_);
  read_attribute(group, "A", A_);
  read_attribute(group, "metastable", metastable_);
  read_attribute(group, "atomic_weight_ratio", awr_);

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
  close_group(kT_group);

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
  double T_min = settings::temperature_range[0];
  double T_max = settings::temperature_range[1];
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

  // Read reactions
  hid_t rxs_group = open_group(group, "reactions");
  for (auto name : group_names(rxs_group)) {
    if (starts_with(name, "reaction_")) {
      hid_t rx_group = open_group(rxs_group, name.c_str());
      reactions_.push_back(std::make_unique<Reaction>(rx_group, temps_to_read));
      close_group(rx_group);
    }
  }
  close_group(rxs_group);

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
  data::nuclides.push_back(std::make_unique<Nuclide>(group, temperature, n));
  return data::nuclides.back().get();
}

extern "C" Reaction* nuclide_reaction(Nuclide* nuc, int i_rx)
{
  return nuc->reactions_[i_rx-1].get();
}

extern "C" void nuclides_clear() { data::nuclides.clear(); }

} // namespace openmc
