#include "openmc/thermal.h"

#include <algorithm> // for sort, move, min, max, find
#include <cmath>     // for round, sqrt, fabs
#include <sstream>   // for stringstream

#include "xtensor/xarray.hpp"
#include "xtensor/xbuilder.hpp"
#include "xtensor/xmath.hpp"
#include "xtensor/xsort.hpp"
#include "xtensor/xtensor.hpp"
#include "xtensor/xview.hpp"

#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/random_lcg.h"
#include "openmc/search.h"
#include "openmc/secondary_correlated.h"
#include "openmc/settings.h"

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
  int secondary_mode;
  if (sec_mode == "equal") {
    secondary_mode = SAB_SECONDARY_EQUAL;
  } else if (sec_mode == "skewed") {
    secondary_mode = SAB_SECONDARY_SKEWED;
  } else if (sec_mode == "continuous") {
    secondary_mode = SAB_SECONDARY_CONT;
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
    data_.emplace_back(T_group, secondary_mode);
    close_group(T_group);
  }

  close_group(kT_group);
}

void
ThermalScattering::calculate_xs(double E, double sqrtkT, int* i_temp,
                                double* elastic, double* inelastic) const
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
  int i_grid = 0;
  double f = 0.0;
  if (E >= sab.inelastic_e_in_.front()) {
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

    if (sab.elastic_mode_ == SAB_ELASTIC_COHERENT) {
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
      f = (E - E_in[i_grid]) / (E_in[i_grid+1] - E_in[i_grid]);

      // Calculate S(a,b) elastic scattering cross section
      auto& xs = sab.elastic_P_;
      *elastic = (1.0 - f) * xs[i_grid] + f * xs[i_grid + 1];
    }
  } else {
    // No elastic data
    *elastic = 0.0;
  }
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

ThermalData::ThermalData(hid_t group, int secondary_mode)
  : inelastic_mode_{secondary_mode}
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

    // Determine whether elastic scattering is incoherent or coherent
    std::string type;
    read_attribute(dset, "type", type);
    if (type == "tab1") {
      elastic_mode_ = SAB_ELASTIC_INCOHERENT;
    } else if (type == "bragg") {
      elastic_mode_ = SAB_ELASTIC_COHERENT;
    }
    close_dataset(dset);

    // Set elastic threshold
    threshold_elastic_ = elastic_e_in_.back();

    // Read angle distribution
    if (elastic_mode_ == SAB_ELASTIC_INCOHERENT) {
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
      n_inelastic_e_out_ = inelastic_e_out_.shape()[1];

      // Read angle distribution
      xt::xarray<double> mu_out;
      read_dataset(inelastic_group, "mu_out", mu_out);
      inelastic_mu_ = mu_out;
      n_inelastic_mu_ = inelastic_mu_.shape()[2];
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

void
ThermalData::sample(const NuclideMicroXS* micro_xs, double E,
                    double* E_out, double* mu)
{
  // Determine whether inelastic or elastic scattering will occur
  if (prn() < micro_xs->thermal_elastic / micro_xs->thermal) {
    // elastic scattering

    // Get index and interpolation factor for elastic grid
    int i = 0;
    double f = 0.0;
    if (E >= elastic_e_in_.front()) {
      auto& E_in = elastic_e_in_;
      i = lower_bound_index(E_in.begin(), E_in.end(), E);
      f = (E - E_in[i]) / (E_in[i+1] - E_in[i]);
    }

    // Select treatment based on elastic mode
    if (elastic_mode_ == SAB_ELASTIC_INCOHERENT) {
      // With this treatment, we interpolate between two discrete cosines
      // corresponding to neighboring incoming energies. This is used for
      // data derived in the incoherent approximation

      // Sample outgoing cosine bin
      int k = prn() * n_elastic_mu_;

      // Determine outgoing cosine corresponding to E_in[i] and E_in[i+1]
      double mu_ijk  = elastic_mu_(i, k);
      double mu_i1jk = elastic_mu_(i+1, k);

      // Cosine of angle between incoming and outgoing neutron
      *mu = (1 - f)*mu_ijk + f*mu_i1jk;

    } else if (elastic_mode_ == SAB_ELASTIC_COHERENT) {
      // This treatment is used for data derived in the coherent
      // approximation, i.e. for crystalline structures that have Bragg
      // edges.

      // Sample a Bragg edge between 1 and i
      double prob = prn() * elastic_P_[i+1];
      int k = 0;
      if (prob >= elastic_P_.front()) {
        k = lower_bound_index(elastic_P_.begin(), elastic_P_.begin() + (i+1), prob);
      }

      // Characteristic scattering cosine for this Bragg edge
      *mu = 1.0 - 2.0*elastic_e_in_[k] / E;

    }

    // Outgoing energy is same as incoming energy
    *E_out = E;

  } else {
    // Perform inelastic calculations

    // Get index and interpolation factor for inelastic grid
    int i = 0;
    double f = 0.0;
    if (E >= inelastic_e_in_.front()) {
      auto& E_in = inelastic_e_in_;
      i = lower_bound_index(E_in.begin(), E_in.end(), E);
      f = (E - E_in[i]) / (E_in[i+1] - E_in[i]);
    }

    // Now that we have an incoming energy bin, we need to determine the
    // outgoing energy bin. This will depend on the "secondary energy
    // mode". If the mode is 0, then the outgoing energy bin is chosen from a
    // set of equally-likely bins. If the mode is 1, then the first
    // two and last two bins are skewed to have lower probabilities than the
    // other bins (0.1 for the first and last bins and 0.4 for the second and
    // second to last bins, relative to a normal bin probability of 1).
    // Finally, if the mode is 2, then a continuous distribution (with
    // accompanying PDF and CDF is utilized)

    if (inelastic_mode_ == SAB_SECONDARY_EQUAL ||
        inelastic_mode_ == SAB_SECONDARY_SKEWED) {
      int j;
      if (inelastic_mode_ == SAB_SECONDARY_EQUAL) {
        // All bins equally likely
        j = prn() * n_inelastic_e_out_;
      } else if (inelastic_mode_ == SAB_SECONDARY_SKEWED) {
        // Distribution skewed away from edge points
        double r = prn() * (n_inelastic_e_out_ - 3);
        if (r > 1.0) {
          // equally likely N-4 middle bins
          j = r + 1;
        } else if (r > 0.6) {
          // second to last bin has relative probability of 0.4
          j = n_inelastic_e_out_ - 2;
        } else if (r > 0.5) {
          // last bin has relative probability of 0.1
          j = n_inelastic_e_out_ - 1;
        } else if (r > 0.1) {
          // second bin has relative probability of 0.4
          j = 1;
        } else {
          // first bin has relative probability of 0.1
          j = 0;
        }
      }

      // Determine outgoing energy corresponding to E_in[i] and E_in[i+1]
      double E_ij  = inelastic_e_out_(i, j);
      double E_i1j = inelastic_e_out_(i+1, j);

      // Outgoing energy
      *E_out = (1 - f)*E_ij + f*E_i1j;

      // Sample outgoing cosine bin
      int k = prn() * n_inelastic_mu_;

      // Determine outgoing cosine corresponding to E_in[i] and E_in[i+1]
      double mu_ijk  = inelastic_mu_(i, j, k);
      double mu_i1jk = inelastic_mu_(i+1, j, k);

      // Cosine of angle between incoming and outgoing neutron
      *mu = (1 - f)*mu_ijk + f*mu_i1jk;

    } else if (inelastic_mode_ == SAB_SECONDARY_CONT) {
      // Continuous secondary energy - this is to be similar to
      // Law 61 interpolation on outgoing energy

      // Sample between ith and [i+1]th bin
      int l = f > prn() ? i + 1 : i;

      // Determine endpoints on grid i
      auto n = inelastic_data_[i].e_out.size();
      double E_i_1 = inelastic_data_[i].e_out(0);
      double E_i_J = inelastic_data_[i].e_out(n);

      // Determine endpoints on grid i + 1
      n = inelastic_data_[i].e_out.size();
      double E_i1_1 = inelastic_data_[i + 1].e_out(0);
      double E_i1_J = inelastic_data_[i + 1].e_out(n);

      double E_1 = E_i_1 + f * (E_i1_1 - E_i_1);
      double E_J = E_i_J + f * (E_i1_J - E_i_J);

      // Determine outgoing energy bin
      // (First reset n_energy_out to the right value)
      n = inelastic_data_[l].n_e_out;
      double r1 = prn();
      double c_j = inelastic_data_[l].e_out_cdf[0];
      double c_j1;
      std::size_t j;
      for (j = 0; j < n - 2; ++j) {
        c_j1 = inelastic_data_[l].e_out_cdf[j + 1];
        if (r1 < c_j1) break;
        c_j = c_j1;
      }

      // check to make sure j is <= n_energy_out - 2
      j = std::min(j, n - 2);

      // Get the data to interpolate between
      double E_l_j = inelastic_data_[l].e_out[j];
      double p_l_j = inelastic_data_[l].e_out_pdf[j];

      // Next part assumes linear-linear interpolation in standard
      double E_l_j1 = inelastic_data_[l].e_out[j + 1];
      double p_l_j1 = inelastic_data_[l].e_out_pdf[j + 1];

      // Find secondary energy (variable E)
      double frac = (p_l_j1 - p_l_j) / (E_l_j1 - E_l_j);
      if (frac == 0.0) {
        *E_out = E_l_j + (r1 - c_j) / p_l_j;
      } else {
        *E_out = E_l_j + (std::sqrt(std::max(0.0, p_l_j*p_l_j +
              2.0*frac*(r1 - c_j))) - p_l_j) / frac;
      }

      // Now interpolate between incident energy bins i and i + 1
      if (l == i) {
        *E_out = E_1 + (E - E_i_1) * (E_J - E_1) / (E_i_J - E_i_1);
      } else {
        *E_out = E_1 + (E - E_i1_1) * (E_J - E_1) / (E_i1_J - E_i1_1);
      }

      // Sample outgoing cosine bin
      std::size_t k = prn() * n_inelastic_mu_;

      // Rather than use the sampled discrete mu directly, it is smeared over
      // a bin of width min(mu[k] - mu[k-1], mu[k+1] - mu[k]) centered on the
      // discrete mu value itself.
      const auto& mu_l = inelastic_data_[l].mu;
      f = (r1 - c_j)/(c_j1 - c_j);

      // Determine (k-1)th mu value
      double mu_left;
      if (k == 0) {
        mu_left = -1.0;
      } else {
        mu_left = mu_l(j, k-1) + f*(mu_l(j+1, k-1) - mu_l(j, k-1));
      }

      // Determine kth mu value
      *mu = mu_l(j, k) + f*(mu_l(j+1, k) - mu_l(j, k));

      // Determine (k+1)th mu value
      double mu_right;
      if (k == n_inelastic_mu_ - 1) {
        mu_right = 1.0 - *mu;
      } else {
        mu_right = mu_l(j, k+1) + f*(mu_l(j+1, k+1) - mu_l(j, k+1)) - *mu;
      }

      // Smear angle
      *mu += std::min(*mu - mu_left, mu_right - *mu)*(prn() - 0.5);

    }  // (inelastic secondary energy treatment)
  }  // (elastic or inelastic)

  // Because of floating-point roundoff, it may be possible for mu to be
  // outside of the range [-1,1). In these cases, we just set mu to exactly
  // -1 or 1
  if (std::fabs(*mu) > 1.0) *mu = std::copysign(1.0, *mu);

}

//==============================================================================
// Fortran compatibility functions
//==============================================================================

ThermalScattering*
sab_from_hdf5(hid_t group, const double* temperature, int n,
    int method, double tolerance, const double* minmax)
{
  // Convert temperatures to a vector
  std::vector<double> T {temperature, temperature + n};

  // Create new object and return it
  return new ThermalScattering{group, T, method, tolerance, minmax};
}

void sab_calculate_xs(ThermalScattering* data, double E, double sqrtkT,
  int* i_temp, double* elastic, double* inelastic)
{
  // Calculate cross section
  int t;
  data->calculate_xs(E, sqrtkT, &t, elastic, inelastic);

  // Fortran needs index plus one
  *i_temp = t + 1;
}

void sab_free(ThermalScattering* data) { delete data; }

bool sab_has_nuclide(ThermalScattering* data, const char* name)
{
  return data->has_nuclide(name);
}

void sab_sample(ThermalScattering* data, const NuclideMicroXS* micro_xs,
  double E_in, double* E_out, double* mu)
{
  int i_temp = micro_xs->index_temp_sab;
  data->data_[i_temp - 1].sample(micro_xs, E_in, E_out, mu);
}

double sab_threshold(ThermalScattering* data) { return data->threshold(); }

} // namespace openmc
