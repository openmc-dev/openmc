#include "thermal.h"

#include <algorithm> // for sort, move
#include <cmath>     // for round
#include <sstream>   // for stringstream

#include "xtensor/xarray.hpp"
#include "xtensor/xbuilder.hpp"
#include "xtensor/xmath.hpp"
#include "xtensor/xsort.hpp"
#include "xtensor/xtensor.hpp"
#include "xtensor/xview.hpp"

#include "constants.h"
#include "error.h"
#include "random_lcg.h"
#include "search.h"
#include "secondary_correlated.h"
#include "settings.h"

namespace openmc {

//==============================================================================
// ThermalScattering implementation
//==============================================================================

ThermalScattering::ThermalScattering(hid_t group, const std::vector<double>& temperature,
                                     int method, double tolerance, const double* minmax)
{
  // Get name of table from group
  name_ = object_name(group);

  // Get rid of leading '/'
  name_ = name_.substr(1);

  read_attribute(group, "atomic_weight_ratio", awr_);
  read_attribute(group, "nuclides", nuclides_);
  std::string sec_mode;
  read_attribute(group, "secondary_mode", sec_mode);
  if (sec_mode == "equal") {
    secondary_mode_ = SAB_SECONDARY_EQUAL;
  } else if (sec_mode == "skewed") {
    secondary_mode_ = SAB_SECONDARY_SKEWED;
  } else if (sec_mode == "continuous") {
    secondary_mode_ = SAB_SECONDARY_CONT;
  }

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
  if (minmax[1] > 0.0) {
    for (const auto& T : temps_available) {
      if (minmax[0] <= T && T <= minmax[1]) {
        temps_to_read.push_back(std::round(T));
      }
    }
  }

  switch (method) {
  case TEMPERATURE_NEAREST:
    // Determine actual temperatures to read
    for (const auto& T : temperature) {

      auto i_closest = xt::argmin(xt::abs(temps_available - T))[0];
      auto temp_actual = temps_available[i_closest];
      if (std::fabs(temp_actual - T) < tolerance) {
        if (std::find(temps_to_read.begin(), temps_to_read.end(), std::round(temp_actual))
            == temps_to_read.end()) {
          temps_to_read.push_back(std::round(temp_actual));
        }
      } else {
        std::stringstream msg;
        msg << "Nuclear data library does not contain cross sections for "
          << name_ << " at or near " << std::round(T) << " K.";
        fatal_error(msg);
      }
    }
    break;

  case TEMPERATURE_INTERPOLATION:
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
        std::stringstream msg;
        msg << "Nuclear data library does not contain cross sections for "
          << name_ << " at temperatures that bound " << std::round(T) << " K.";
        fatal_error(msg);
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
    std::string temp_str = std::to_string(T) + "K";

    // Read exact temperature value
    double kT;
    read_dataset(kT_group, temp_str.data(), kT);
    kTs_.push_back(kT);

    // Open group for temperature i
    hid_t T_group = open_group(group, temp_str.data());
    data_.emplace_back(T_group, secondary_mode_);
    close_group(group);
  }

  close_group(kT_group);
}

void
ThermalScattering::calculate_xs(double E, double sqrtkT, int* i_temp,
                                double* elastic, double* inelastic)
{
  // Determine temperature for S(a,b) table
  double kT = sqrtkT*sqrtkT;
  int i;
  if (temperature_method == TEMPERATURE_NEAREST) {
    // If using nearest temperature, do linear search on temperature
    for (i = 0; i < kTs_.size(); ++i) {
      if (abs(kTs_[i] - kT) < K_BOLTZMANN*temperature_tolerance) {
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
    if (f > prn()) ++i;
  }

  // Set temperature index
  *i_temp = i;

  // Get pointer to S(a,b) table
  auto& sab = data_[i];

  // Get index and interpolation factor for inelastic grid
  int i_grid;
  double f;
  if (E < sab.inelastic_e_in_.front()) {
    i_grid = 0;
    f = 0.0;
  } else {
    auto& E_in = sab.inelastic_e_in_;
    i_grid = lower_bound_index(E_in.begin(), E_in.end(), E);
    f = (E - E_in[i_grid]) / (E_in[i_grid+1] - E_in[i_grid]);
  }

  // Calculate S(a,b) inelastic scattering cross section
  auto& xs = sab.inelastic_sigma_;
  *inelastic = (1.0 - f) * xs[i_grid] + f * xs[i_grid + 1];

  // Check for elastic data
  if (E < sab.threshold_elastic_) {
    // Determine whether elastic scattering is given in the coherent or
    // incoherent approximation. For coherent, the cross section is
    // represented as P/E whereas for incoherent, it is simply P

    auto& E_in = sab.elastic_e_in_;

    if (sab.elastic_mode_ == SAB_ELASTIC_EXACT) {
      if (E < E_in.front()) {
        // If energy is below that of the lowest Bragg peak, the elastic
        // cross section will be zero
        *elastic = 0.0;
      } else {
        i_grid = lower_bound_index(E_in.begin(), E_in.end(), E);
        *elastic = sab.elastic_P_[i_grid] / E;
      }
    } else {
      // Determine index on elastic energy grid
      if (E < E_in.front()) {
        i_grid = 0;
      } else {
        i_grid = lower_bound_index(E_in.begin(), E_in.end(), E);
      }

      // Get interpolation factor for elastic grid
      f = (E - E_in[i_grid])/(E_in[i_grid+1] - E_in[i_grid]);

      // Calculate S(a,b) elastic scattering cross section
      auto& xs = sab.elastic_P_;
      *elastic = (1.0 - f) * xs[i_grid] + f * xs[i_grid + 1];
    }
  } else {
    // No elastic data
    *elastic = 0.0;
  }
}

//==============================================================================
// ThermalData implementation
//==============================================================================

ThermalData::ThermalData(hid_t group, int secondary_mode)
{
  // Coherent elastic data
  if (object_exists(group, "elastic")) {
    // Read cross section data
    hid_t elastic_group = open_group(group, "elastic");

    // Read elastic cross section
    xt::xarray<double> temp;
    hid_t dset = open_dataset(elastic_group, "xs");
    read_dataset(dset, temp);

    // Get view on energies and cross section/probability values
    auto E_in = xt::view(temp, 0);
    auto P = xt::view(temp, 1);

    // Set cross section data and type
    std::copy(E_in.begin(), E_in.end(), std::back_inserter(elastic_e_in_));
    std::copy(P.begin(), P.end(), std::back_inserter(elastic_P_));
    n_elastic_e_in_ = elastic_e_in_.size();

    // Determine elastic type
    std::string type;
    read_attribute(dset, "type", type);
    if (type == "tab1") {
      elastic_mode_ = SAB_ELASTIC_DISCRETE;
    } else if (type == "bragg") {
      elastic_mode_ = SAB_ELASTIC_EXACT;
    }
    close_dataset(dset);

    // Set elastic threshold
    threshold_elastic_ = elastic_e_in_.back();

    // Read angle distribution
    if (elastic_mode_ != SAB_ELASTIC_EXACT) {
      xt::xarray<double> mu_out;
      read_dataset(elastic_group, "mu_out", mu_out);
      elastic_mu_ = mu_out;
    }

    close_group(elastic_group);
  }

  // Inelastic data
  if (object_exists(group, "inelastic")) {
    // Read type of inelastic data
    hid_t inelastic_group = open_group(group, "inelastic");

    // Read cross section data
    xt::xarray<double> temp;
    read_dataset(inelastic_group, "xs", temp);

    // Get view of inelastic cross section and energy grid
    auto E_in = xt::view(temp, 0);
    auto xs = xt::view(temp, 1);

    // Set cross section data
    std::copy(E_in.begin(), E_in.end(), std::back_inserter(inelastic_e_in_));
    std::copy(xs.begin(), xs.end(), std::back_inserter(inelastic_sigma_));
    n_inelastic_e_in_ = inelastic_e_in_.size();

    // Set inelastic threshold
    threshold_inelastic_ = inelastic_e_in_.back();

    if (secondary_mode != SAB_SECONDARY_CONT) {
      // Read energy distribution
      xt::xarray<double> E_out;
      read_dataset(inelastic_group, "energy_out", E_out);
      inelastic_e_out_ = E_out;

      // Read angle distribution
      xt::xarray<double> mu_out;
      read_dataset(inelastic_group, "mu_out", mu_out);
      inelastic_mu_ = mu_out;
    } else {
      // Read correlated angle-energy distribution
      CorrelatedAngleEnergy dist {inelastic_group};

      // Convert to S(a,b) native format
      for (const auto& edist : dist.distribution()) {
        // Create temporary distribution
        DistEnergySab d;

        // Copy outgoing energy distribution
        d.n_e_out = edist.e_out.size();
        d.e_out = edist.e_out;
        d.e_out_pdf = edist.p;
        d.e_out_cdf = edist.c;

        for (int j = 0; j < d.n_e_out; ++j) {
          auto adist = dynamic_cast<Tabular*>(edist.angle[j].get());
          if (adist) {
            // On first pass, allocate space for angles
            if (j == 0) {
              auto n_mu = adist->x().size();
              n_inelastic_mu_ = n_mu;
              d.mu = xt::empty<double>({d.n_e_out, n_mu});
            }

            // Copy outgoing angles
            auto mu_j = xt::view(d.mu, j);
            std::copy(adist->x().begin(), adist->x().end(), mu_j.begin());
          }
        }

        inelastic_data_.push_back(std::move(d));
      }
    }

    close_group(inelastic_group);
  }
}

} // namespace openmc
