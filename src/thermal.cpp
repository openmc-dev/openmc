#include "openmc/thermal.h"

#include <algorithm> // for sort, move, min, max, find
#include <cmath>     // for round, sqrt, abs

#include <fmt/core.h>
#include "xtensor/xarray.hpp"
#include "xtensor/xbuilder.hpp"
#include "xtensor/xmath.hpp"
#include "xtensor/xsort.hpp"
#include "xtensor/xtensor.hpp"
#include "xtensor/xview.hpp"

#include "openmc/constants.h"
#include "openmc/endf.h"
#include "openmc/error.h"
#include "openmc/random_lcg.h"
#include "openmc/search.h"
#include "openmc/secondary_correlated.h"
#include "openmc/secondary_thermal.h"
#include "openmc/settings.h"

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace data {
std::vector<std::unique_ptr<ThermalScattering>> thermal_scatt;
std::unordered_map<std::string, int> thermal_scatt_map;
}

//==============================================================================
// ThermalScattering implementation
//==============================================================================

ThermalScattering::ThermalScattering(hid_t group, const std::vector<double>& temperature)
{
  // Get name of table from group
  name_ = object_name(group);

  // Get rid of leading '/'
  name_ = name_.substr(1);

  read_attribute(group, "atomic_weight_ratio", awr_);
  read_attribute(group, "energy_max", energy_max_);
  read_attribute(group, "nuclides", nuclides_);

  // Read temperatures
  hid_t kT_group = open_group(group, "kTs");

  // Determine temperatures available
  auto dset_names = dataset_names(kT_group);
  auto n = dset_names.size();
  auto temps_available = xt::empty<double>({n});
  for (int i = 0; i < dset_names.size(); ++i) {
    // Read temperature value
    double T;
    read_dataset(kT_group, dset_names[i].data(), T);
    temps_available[i] = T / K_BOLTZMANN;
  }
  std::sort(temps_available.begin(), temps_available.end());

  // Determine actual temperatures to read -- start by checking whether a
  // temperature range was given, in which case all temperatures in the range
  // are loaded irrespective of what temperatures actually appear in the model
  std::vector<int> temps_to_read;
  if (settings::temperature_range[1] > 0.0) {
    for (const auto& T : temps_available) {
      if (settings::temperature_range[0] <= T &&
          T <= settings::temperature_range[1]) {
        temps_to_read.push_back(std::round(T));
      }
    }
  }

  switch (settings::temperature_method) {
  case TemperatureMethod::NEAREST:
    // Determine actual temperatures to read
    for (const auto& T : temperature) {

      auto i_closest = xt::argmin(xt::abs(temps_available - T))[0];
      auto temp_actual = temps_available[i_closest];
      if (std::abs(temp_actual - T) < settings::temperature_tolerance) {
        if (std::find(temps_to_read.begin(), temps_to_read.end(), std::round(temp_actual))
            == temps_to_read.end()) {
          temps_to_read.push_back(std::round(temp_actual));
        }
      } else {
        fatal_error(fmt::format("Nuclear data library does not contain cross "
          "sections for {} at or near {} K.", name_, std::round(T)));
      }
    }
    break;

  case TemperatureMethod::INTERPOLATION:
    // If temperature interpolation or multipole is selected, get a list of
    // bounding temperatures for each actual temperature present in the model
    for (const auto& T : temperature) {
      bool found = false;
      for (int j = 0; j < temps_available.size() - 1; ++j) {
        if (temps_available[j] <= T && T < temps_available[j + 1]) {
          int T_j = std::round(temps_available[j]);
          int T_j1 = std::round(temps_available[j + 1]);
          if (std::find(temps_to_read.begin(), temps_to_read.end(), T_j) == temps_to_read.end()) {
            temps_to_read.push_back(T_j);
          }
          if (std::find(temps_to_read.begin(), temps_to_read.end(), T_j1) == temps_to_read.end()) {
            temps_to_read.push_back(T_j1);
          }
          found = true;
        }
      }
      if (!found) {
        fatal_error(fmt::format("Nuclear data library does not contain cross "
          "sections for {} at temperatures that bound {} K.", name_, std::round(T)));
      }
    }
  }

  // Sort temperatures to read
  std::sort(temps_to_read.begin(), temps_to_read.end());

  auto n_temperature = temps_to_read.size();
  kTs_.reserve(n_temperature);
  data_.reserve(n_temperature);

  for (auto T : temps_to_read) {
    // Get temperature as a string
    std::string temp_str = fmt::format("{}K", T);

    // Read exact temperature value
    double kT;
    read_dataset(kT_group, temp_str.data(), kT);
    kTs_.push_back(kT);

    // Open group for temperature i
    hid_t T_group = open_group(group, temp_str.data());
    data_.emplace_back(T_group);
    close_group(T_group);
  }

  close_group(kT_group);
}

void
ThermalScattering::calculate_xs(double E, double sqrtkT, int* i_temp,
                                double* elastic, double* inelastic,
                                uint64_t* seed) const
{
  // Determine temperature for S(a,b) table
  double kT = sqrtkT*sqrtkT;
  int i;
  if (settings::temperature_method == TemperatureMethod::NEAREST) {
    // If using nearest temperature, do linear search on temperature
    for (i = 0; i < kTs_.size(); ++i) {
      if (std::abs(kTs_[i] - kT) < K_BOLTZMANN*settings::temperature_tolerance) {
        break;
      }
    }
  } else {
    // Find temperatures that bound the actual temperature
    for (i = 0; i < kTs_.size() - 1; ++i) {
      if (kTs_[i] <= kT && kT < kTs_[i+1]) {
        break;
      }
    }

    // Randomly sample between temperature i and i+1
    double f = (kT - kTs_[i]) / (kTs_[i+1] - kTs_[i]);
    if (f > prn(seed)) ++i;
  }

  // Set temperature index
  *i_temp = i;

  // Calculate cross sections for ith temperature
  data_[i].calculate_xs(E, elastic, inelastic);
}

bool
ThermalScattering::has_nuclide(const char* name) const
{
  std::string nuc {name};
  return std::find(nuclides_.begin(), nuclides_.end(), nuc) != nuclides_.end();
}

//==============================================================================
// ThermalData implementation
//==============================================================================

ThermalData::ThermalData(hid_t group)
{
  // Coherent/incoherent elastic data
  if (object_exists(group, "elastic")) {
    // Read cross section data
    hid_t elastic_group = open_group(group, "elastic");

    // Read elastic cross section
    elastic_.xs = read_function(elastic_group, "xs");

    // Read angle-energy distribution
    hid_t dgroup = open_group(elastic_group, "distribution");
    std::string temp;
    read_attribute(dgroup, "type", temp);
    if (temp == "coherent_elastic") {
      auto xs = dynamic_cast<CoherentElasticXS*>(elastic_.xs.get());
      elastic_.distribution = std::make_unique<CoherentElasticAE>(*xs);
    } else {
      if (temp == "incoherent_elastic") {
        elastic_.distribution = std::make_unique<IncoherentElasticAE>(dgroup);
      } else if (temp == "incoherent_elastic_discrete") {
        auto xs = dynamic_cast<Tabulated1D*>(elastic_.xs.get());
        elastic_.distribution = std::make_unique<IncoherentElasticAEDiscrete>(
          dgroup, xs->x()
        );
      }
    }

    close_group(elastic_group);
  }

  // Inelastic data
  if (object_exists(group, "inelastic")) {
    // Read type of inelastic data
    hid_t inelastic_group = open_group(group, "inelastic");

    // Read inelastic cross section
    inelastic_.xs = read_function(inelastic_group, "xs");

    // Read angle-energy distribution
    hid_t dgroup = open_group(inelastic_group, "distribution");
    std::string temp;
    read_attribute(dgroup, "type", temp);
    if (temp == "incoherent_inelastic") {
      inelastic_.distribution = std::make_unique<IncoherentInelasticAE>(dgroup);
    } else if (temp == "incoherent_inelastic_discrete") {
      auto xs = dynamic_cast<Tabulated1D*>(inelastic_.xs.get());
      inelastic_.distribution = std::make_unique<IncoherentInelasticAEDiscrete>(
        dgroup, xs->x()
      );
    }

    close_group(inelastic_group);
  }
}

void
ThermalData::calculate_xs(double E, double* elastic, double* inelastic) const
{
  // Calculate thermal elastic scattering cross section
  if (elastic_.xs) {
    *elastic = (*elastic_.xs)(E);
  } else {
    *elastic = 0.0;
  }

  // Calculate thermal inelastic scattering cross section
  *inelastic = (*inelastic_.xs)(E);
}

void
ThermalData::sample(const NuclideMicroXS& micro_xs, double E,
                    double* E_out, double* mu, uint64_t* seed)
{
  // Determine whether inelastic or elastic scattering will occur
  if (prn(seed) < micro_xs.thermal_elastic / micro_xs.thermal) {
    elastic_.distribution->sample(E, *E_out, *mu, seed);
  } else {
    inelastic_.distribution->sample(E, *E_out, *mu, seed);
  }

  // Because of floating-point roundoff, it may be possible for mu to be
  // outside of the range [-1,1). In these cases, we just set mu to exactly
  // -1 or 1
  if (std::abs(*mu) > 1.0) *mu = std::copysign(1.0, *mu);
}

void free_memory_thermal()
{
  data::thermal_scatt.clear();
  data::thermal_scatt_map.clear();
}

} // namespace openmc
