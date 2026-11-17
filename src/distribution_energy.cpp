#include "openmc/distribution_energy.h"

#include <algorithm> // for max, min, copy, move
#include <cmath>     // for isfinite
#include <cstddef>   // for size_t
#include <iterator>  // for back_inserter

#include "openmc/tensor.h"

#include "openmc/endf.h"
#include "openmc/hdf5_interface.h"
#include "openmc/math_functions.h"
#include "openmc/random_dist.h"
#include "openmc/random_lcg.h"
#include "openmc/search.h"

namespace openmc {

//==============================================================================
// DiscretePhoton implementation
//==============================================================================

DiscretePhoton::DiscretePhoton(hid_t group)
{
  read_attribute(group, "primary_flag", primary_flag_);
  read_attribute(group, "energy", energy_);
  read_attribute(group, "atomic_weight_ratio", A_);
}

double DiscretePhoton::sample(double E, uint64_t* seed) const
{
  if (primary_flag_ == 2) {
    return energy_ + A_ / (A_ + 1) * E;
  } else {
    return energy_;
  }
}

//==============================================================================
// LevelInelastic implementation
//==============================================================================

LevelInelastic::LevelInelastic(hid_t group)
{
  read_attribute(group, "threshold", threshold_);
  read_attribute(group, "mass_ratio", mass_ratio_);
}

double LevelInelastic::sample(double E, uint64_t* seed) const
{
  return mass_ratio_ * (E - threshold_);
}

//==============================================================================
// ContinuousTabular implementation
//==============================================================================

ContinuousTabular::ContinuousTabular(hid_t group)
{
  // Open incoming energy dataset
  hid_t dset = open_dataset(group, "energy");

  // Get interpolation parameters
  tensor::Tensor<int> temp;
  read_attribute(dset, "interpolation", temp);

  tensor::View<int> temp_b = temp.slice(0); // breakpoints
  tensor::View<int> temp_i = temp.slice(1); // interpolation parameters

  std::copy(temp_b.begin(), temp_b.end(), std::back_inserter(breakpoints_));
  for (const auto i : temp_i)
    interpolation_.push_back(int2interp(i));
  n_region_ = breakpoints_.size();

  // Get incoming energies
  read_dataset(dset, energy_);
  std::size_t n_energy = energy_.size();
  close_dataset(dset);

  // Get outgoing energy distribution data
  dset = open_dataset(group, "distribution");
  vector<int> offsets;
  vector<int> interp;
  vector<int> n_discrete;
  read_attribute(dset, "offsets", offsets);
  read_attribute(dset, "interpolation", interp);
  read_attribute(dset, "n_discrete_lines", n_discrete);

  tensor::Tensor<double> eout;
  read_dataset(dset, eout);
  close_dataset(dset);

  for (int i = 0; i < n_energy; ++i) {
    // Determine number of outgoing energies
    int j = offsets[i];
    int n;
    if (i < n_energy - 1) {
      n = offsets[i + 1] - j;
    } else {
      n = eout.shape(1) - j;
    }

    // Assign interpolation scheme and number of discrete lines
    CTTable d;
    d.interpolation = int2interp(interp[i]);
    d.n_discrete = n_discrete[i];

    // Copy data
    d.e_out = eout.slice(0, tensor::range(j, j + n));
    d.p = eout.slice(1, tensor::range(j, j + n));

    // To get answers that match ACE data, for now we still use the tabulated
    // CDF values that were passed through to the HDF5 library. At a later
    // time, we can remove the CDF values from the HDF5 library and
    // reconstruct them using the PDF
    if (true) {
      d.c = eout.slice(2, tensor::range(j, j + n));
    } else {
      // Calculate cumulative distribution function -- discrete portion
      for (int k = 0; k < d.n_discrete; ++k) {
        if (k == 0) {
          d.c[k] = d.p[k];
        } else {
          d.c[k] = d.c[k - 1] + d.p[k];
        }
      }

      // Continuous portion
      for (int k = d.n_discrete; k < n; ++k) {
        if (k == d.n_discrete) {
          d.c[k] = d.c[k - 1] + d.p[k];
        } else {
          if (d.interpolation == Interpolation::histogram) {
            d.c[k] = d.c[k - 1] + d.p[k - 1] * (d.e_out[k] - d.e_out[k - 1]);
          } else if (d.interpolation == Interpolation::lin_lin) {
            d.c[k] = d.c[k - 1] + 0.5 * (d.p[k - 1] + d.p[k]) *
                                    (d.e_out[k] - d.e_out[k - 1]);
          }
        }
      }

      // Normalize density and distribution functions
      d.p /= d.c[n - 1];
      d.c /= d.c[n - 1];
    }

    distribution_.push_back(std::move(d));
  } // incoming energies
}

void ContinuousTabular::set_spectrum(const vector<double>& energy,
  const vector<int>& interpolation, const vector<int>& n_discrete,
  const vector<vector<double>>& e_out, const vector<vector<double>>& p)
{
  const auto n = energy_.size();
  if (energy != energy_ || n < 2 || distribution_.size() != n ||
      interpolation.size() != n || n_discrete.size() != n ||
      e_out.size() != n || p.size() != n) {
    throw std::invalid_argument("Incident-energy grid and table count must "
                                "exactly match the native spectrum.");
  }
  vector<CTTable> tables;
  tables.reserve(n);
  for (size_t i = 0; i < n; ++i) {
    const auto& original = distribution_[i];
    if (n_discrete[i] == original.n_discrete &&
        interpolation[i] == static_cast<int>(original.interpolation) &&
        e_out[i] ==
          vector<double>(original.e_out.begin(), original.e_out.end()) &&
        p[i] == vector<double>(original.p.begin(), original.p.end())) {
      // Unmodified tables retain their loaded CDF exactly in partial edits.
      tables.push_back(original);
      continue;
    }
    tables.push_back(make_table(static_cast<Interpolation>(interpolation[i]),
      n_discrete[i], e_out[i], p[i],
      n_discrete[i] == distribution_[i].n_discrete &&
        interpolation[i] == static_cast<int>(distribution_[i].interpolation) &&
        e_out[i] == vector<double>(distribution_[i].e_out.begin(),
                      distribution_[i].e_out.end())));
  }

  // Unit-base interpolation needs continuous endpoints in both adjacent tables.
  // Histogram incident interpolation only samples one table and needs no
  // mapping.
  const bool histogram =
    n_region_ == 1 && interpolation_[0] == Interpolation::histogram;
  if (!histogram) {
    for (size_t i = 1; i < n; ++i) {
      bool discrete = tables[i].n_discrete == tables[i].e_out.size();
      bool previous_discrete =
        tables[i - 1].n_discrete == tables[i - 1].e_out.size();
      if (discrete != previous_discrete) {
        throw std::invalid_argument("Unit-base interpolation cannot mix "
                                    "discrete-only and continuous tables.");
      }
    }
  }
  distribution_.swap(tables);
}

ContinuousTabular::CTTable ContinuousTabular::make_table(
  Interpolation interpolation, int n_discrete, const vector<double>& e_out,
  const vector<double>& p, bool native_grid)
{
  if (e_out.size() != p.size() || e_out.empty()) {
    throw std::invalid_argument("Outgoing-energy and probability arrays must "
                                "be nonempty and equal in length.");
  }
  if (n_discrete < 0 || n_discrete > e_out.size()) {
    throw std::invalid_argument("Invalid number of discrete lines.");
  }
  if (n_discrete < e_out.size() && e_out.size() - n_discrete < 2) {
    throw std::invalid_argument(
      "A continuous outgoing-energy distribution needs at least two points.");
  }
  if (interpolation != Interpolation::histogram &&
      interpolation != Interpolation::lin_lin) {
    throw std::invalid_argument(
      "Only histogram and linear-linear interpolation are supported.");
  }

  double norm = 0.0;
  for (int k = 0; k < e_out.size(); ++k) {
    if (!std::isfinite(e_out[k]) || e_out[k] < 0.0 || !std::isfinite(p[k]) ||
        p[k] < 0.0) {
      throw std::invalid_argument(
        "Outgoing energies and probabilities must be finite and nonnegative.");
    }
    if (k < n_discrete) {
      norm += p[k];
    } else if (k > n_discrete) {
      double width = e_out[k] - e_out[k - 1];
      if (width < 0.0 || (width == 0.0 && !native_grid)) {
        throw std::invalid_argument(
          "Continuous outgoing energies must be strictly increasing.");
      }
      norm += interpolation == Interpolation::histogram
                ? p[k - 1] * width
                : 0.5 * (p[k - 1] + p[k]) * width;
    }
  }
  if (!(norm > 0.0) || !std::isfinite(norm)) {
    throw std::invalid_argument(
      "Outgoing-energy distribution has zero or invalid normalization.");
  }

  CTTable table;
  table.interpolation = interpolation;
  table.n_discrete = n_discrete;
  table.e_out = tensor::Tensor<double>(e_out.data(), e_out.size());
  table.p = tensor::Tensor<double>(p.data(), p.size());
  table.c = tensor::zeros<double>({p.size()});

  double c = 0.0;
  for (int k = 0; k < p.size(); ++k) {
    table.p[k] /= norm;
    if (!std::isfinite(table.p[k])) {
      throw std::invalid_argument("Normalized probability is not finite.");
    }
    if (k < n_discrete) {
      c += table.p[k];
    } else if (k > n_discrete) {
      double width = table.e_out[k] - table.e_out[k - 1];
      c += interpolation == Interpolation::histogram
             ? table.p[k - 1] * width
             : 0.5 * (table.p[k - 1] + table.p[k]) * width;
    }
    table.c[k] = std::min(c, 1.0);
  }
  table.c[table.c.size() - 1] = 1.0;
  return table;
}

double ContinuousTabular::sample(double E, uint64_t* seed) const
{
  // Read number of interpolation regions and incoming energies
  bool histogram_interp;
  if (n_region_ == 1) {
    histogram_interp = (interpolation_[0] == Interpolation::histogram);
  } else {
    histogram_interp = false;
  }

  // Find energy bin and calculate interpolation factor -- if the energy is
  // outside the range of the tabulated energies, choose the first or last bins
  auto n_energy_in = energy_.size();
  int i;
  double r;
  if (E < energy_[0]) {
    i = 0;
    r = 0.0;
  } else if (E > energy_[n_energy_in - 1]) {
    i = n_energy_in - 2;
    r = 1.0;
  } else {
    i = lower_bound_index(energy_.begin(), energy_.end(), E);
    r = (E - energy_[i]) / (energy_[i + 1] - energy_[i]);
  }

  // Sample between the ith and [i+1]th bin
  int l;
  if (histogram_interp) {
    l = i;
  } else {
    l = r > prn(seed) ? i + 1 : i;
  }

  // Determine outgoing energy bin
  int n_energy_out = distribution_[l].e_out.size();
  int n_discrete = distribution_[l].n_discrete;
  double r1 = prn(seed);
  double c_k = distribution_[l].c[0];
  int k = 0;
  int end = n_energy_out - 2;

  // A discrete line is sampled directly, irrespective of the interpolation
  // assigned to the continuous part. In particular, do not access k + 1 for
  // a discrete-only table.
  for (int j = 0; j < n_discrete; ++j) {
    if (r1 < distribution_[l].c[j]) {
      return distribution_[l].e_out[j];
    }
  }
  if (n_discrete == n_energy_out) {
    return distribution_[l].e_out[n_energy_out - 1];
  }
  k = n_discrete;
  c_k = distribution_[l].c[k];

  // Continuous portion
  double c_k1;
  for (int j = n_discrete; j < end; ++j) {
    k = j;
    c_k1 = distribution_[l].c[k + 1];
    if (r1 < c_k1)
      break;
    k = j + 1;
    c_k = c_k1;
  }

  double E_l_k = distribution_[l].e_out[k];
  double p_l_k = distribution_[l].p[k];
  double E_out = E_l_k;
  if (distribution_[l].interpolation == Interpolation::histogram) {
    // Histogram interpolation
    if (p_l_k > 0.0 && k >= n_discrete) {
      E_out = E_l_k + (r1 - c_k) / p_l_k;
    }

  } else if (distribution_[l].interpolation == Interpolation::lin_lin) {
    // Linear-linear interpolation
    double E_l_k1 = distribution_[l].e_out[k + 1];
    double p_l_k1 = distribution_[l].p[k + 1];

    if (E_l_k != E_l_k1) {
      double frac = (p_l_k1 - p_l_k) / (E_l_k1 - E_l_k);
      if (frac == 0.0) {
        E_out = E_l_k + (r1 - c_k) / p_l_k;
      } else {
        E_out =
          E_l_k +
          (std::sqrt(std::max(0.0, p_l_k * p_l_k + 2.0 * frac * (r1 - c_k))) -
            p_l_k) /
            frac;
      }
    }
  } else {
    throw std::runtime_error {"Unexpected interpolation for continuous energy "
                              "distribution."};
  }

  // Now interpolate between incident energy bins i and i + 1
  if (!histogram_interp && n_energy_out > 1 && k >= n_discrete) {
    // Interpolation for energy E1 and EK
    n_energy_out = distribution_[i].e_out.size();
    n_discrete = distribution_[i].n_discrete;
    const double E_i_1 = distribution_[i].e_out[n_discrete];
    const double E_i_K = distribution_[i].e_out[n_energy_out - 1];

    n_energy_out = distribution_[i + 1].e_out.size();
    n_discrete = distribution_[i + 1].n_discrete;
    const double E_i1_1 = distribution_[i + 1].e_out[n_discrete];
    const double E_i1_K = distribution_[i + 1].e_out[n_energy_out - 1];

    const double E_1 = E_i_1 + r * (E_i1_1 - E_i_1);
    const double E_K = E_i_K + r * (E_i1_K - E_i_K);

    if (l == i) {
      return E_1 + (E_out - E_i_1) * (E_K - E_1) / (E_i_K - E_i_1);
    } else {
      return E_1 + (E_out - E_i1_1) * (E_K - E_1) / (E_i1_K - E_i1_1);
    }
  } else {
    return E_out;
  }
}

//==============================================================================
// MaxwellEnergy implementation
//==============================================================================

MaxwellEnergy::MaxwellEnergy(hid_t group)
{
  read_attribute(group, "u", u_);
  hid_t dset = open_dataset(group, "theta");
  theta_ = Tabulated1D {dset};
  close_dataset(dset);
}

double MaxwellEnergy::sample(double E, uint64_t* seed) const
{
  // Get temperature corresponding to incoming energy
  double theta = theta_(E);

  while (true) {
    // Sample maxwell fission spectrum
    double E_out = maxwell_spectrum(theta, seed);

    // Accept energy based on restriction energy
    if (E_out <= E - u_)
      return E_out;
  }
}

//==============================================================================
// Evaporation implementation
//==============================================================================

Evaporation::Evaporation(hid_t group)
{
  read_attribute(group, "u", u_);
  hid_t dset = open_dataset(group, "theta");
  theta_ = Tabulated1D {dset};
  close_dataset(dset);
}

double Evaporation::sample(double E, uint64_t* seed) const
{
  // Get temperature corresponding to incoming energy
  double theta = theta_(E);

  double y = (E - u_) / theta;
  double v = 1.0 - std::exp(-y);

  // Sample outgoing energy based on evaporation spectrum probability
  // density function
  double x;
  while (true) {
    x = -std::log((1.0 - v * prn(seed)) * (1.0 - v * prn(seed)));
    if (x <= y)
      break;
  }

  return x * theta;
}

//==============================================================================
// WattEnergy implementation
//==============================================================================

WattEnergy::WattEnergy(hid_t group)
{
  // Read restriction energy
  read_attribute(group, "u", u_);

  // Read tabulated functions
  hid_t dset = open_dataset(group, "a");
  a_ = Tabulated1D {dset};
  close_dataset(dset);
  dset = open_dataset(group, "b");
  b_ = Tabulated1D {dset};
  close_dataset(dset);
}

double WattEnergy::sample(double E, uint64_t* seed) const
{
  // Determine Watt parameters at incident energy
  double a = a_(E);
  double b = b_(E);

  while (true) {
    // Sample energy-dependent Watt fission spectrum
    double E_out = watt_spectrum(a, b, seed);

    // Accept energy based on restriction energy
    if (E_out <= E - u_)
      return E_out;
  }
}

} // namespace openmc
