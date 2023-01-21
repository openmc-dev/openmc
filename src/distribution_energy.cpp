#include "openmc/distribution_energy.h"

#include <algorithm> // for max, min, copy, move
#include <cstddef>   // for size_t
#include <iterator>  // for back_inserter

#include <iostream>
#include "xtensor/xview.hpp"

#include "openmc/endf.h"
#include "openmc/error.h"
#include "openmc/hdf5_interface.h"
#include "openmc/math_functions.h"
#include "openmc/random_lcg.h"
#include "openmc/search.h"

namespace openmc {

EnergyDistributionFlat::EnergyDistributionFlat(const uint8_t* data)
  : data_{data}
{
  type_ = static_cast<EnergyDistType>(*reinterpret_cast<const int*>(data_));
}

double EnergyDistributionFlat::sample(double E, uint64_t* seed) const
{
  switch (type_) {
  case EnergyDistType::DISCRETE_PHOTON:
    {
      DiscretePhotonFlat dist(data_);
      return dist.sample(E, seed);
    }
  case EnergyDistType::LEVEL_INELASTIC:
    {
      LevelInelasticFlat dist(data_);
      return dist.sample(E, seed);
    }
  case EnergyDistType::CONTINUOUS_TABULAR:
    {
      ContinuousTabularFlat dist(data_);
      return dist.sample(E, seed);
    }
  case EnergyDistType::EVAPORATION:
    {
      EvaporationFlat dist(data_);
      return dist.sample(E, seed);
    }
  case EnergyDistType::MAXWELL:
    {
      MaxwellFlat dist(data_);
      return dist.sample(E, seed);
    }
  case EnergyDistType::WATT:
    {
      WattFlat dist(data_);
      return dist.sample(E, seed);
    }
  default:
    UNREACHABLE();
  }
}

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
    return energy_ + A_/(A_+ 1)*E;
  } else {
    return energy_;
  }
}

void DiscretePhoton::serialize(DataBuffer& buffer) const
{
  buffer.add(static_cast<int>(EnergyDistType::DISCRETE_PHOTON)); // 4
  buffer.add(primary_flag_);                                     // 4
  buffer.add(energy_);                                           // 8
  buffer.add(A_);                                                // 8
}

double DiscretePhotonFlat::sample(double E, uint64_t* seed) const
{
  if (this->primary_flag() == 2) {
    double A = this->A();
    return energy() + A/(A+ 1)*E;
  } else {
    return energy();
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
  return mass_ratio_*(E - threshold_);
}

void LevelInelastic::serialize(DataBuffer& buffer) const
{
  buffer.add(static_cast<int>(EnergyDistType::LEVEL_INELASTIC)); // 4
  buffer.align(8);                                               // 4
  buffer.add(threshold_);                                        // 8
  buffer.add(mass_ratio_);                                       // 8
}

double LevelInelasticFlat::sample(double E, uint64_t* seed) const
{
  return this->mass_ratio()*(E - this->threshold());
}

//==============================================================================
// ContinuousTabular implementation
//==============================================================================

ContinuousTabular::ContinuousTabular(hid_t group)
{
  // Open incoming energy dataset
  hid_t dset = open_dataset(group, "energy");

  // Get interpolation parameters
  xt::xarray<int> temp;
  read_attribute(dset, "interpolation", temp);

  auto temp_b = xt::view(temp, 0); // view of breakpoints
  auto temp_i = xt::view(temp, 1); // view of interpolation parameters

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
  std::vector<int> offsets;
  std::vector<int> interp;
  std::vector<int> n_discrete;
  read_attribute(dset, "offsets", offsets);
  read_attribute(dset, "interpolation", interp);
  read_attribute(dset, "n_discrete_lines", n_discrete);

  xt::xarray<double> eout;
  read_dataset(dset, eout);
  close_dataset(dset);

  for (int i = 0; i < n_energy; ++i) {
    // Determine number of outgoing energies
    int j = offsets[i];
    int n;
    if (i < n_energy - 1) {
      n = offsets[i+1] - j;
    } else {
      n = eout.shape()[1] - j;
    }

    // Assign interpolation scheme and number of discrete lines
    CTTable d;
    d.interpolation = int2interp(interp[i]);
    d.n_discrete = n_discrete[i];

    // Copy data
    d.e_out = xt::view(eout, 0, xt::range(j, j+n));
    d.p = xt::view(eout, 1, xt::range(j, j+n));

    // To get answers that match ACE data, for now we still use the tabulated
    // CDF values that were passed through to the HDF5 library. At a later
    // time, we can remove the CDF values from the HDF5 library and
    // reconstruct them using the PDF
    if (true) {
      d.c = xt::view(eout, 2, xt::range(j, j+n));
    } else {
      // Calculate cumulative distribution function -- discrete portion
      for (int k = 0; k < d.n_discrete; ++k) {
        if (k == 0) {
          d.c[k] = d.p[k];
        } else {
          d.c[k] = d.c[k-1] + d.p[k];
        }
      }

      // Continuous portion
      for (int k = d.n_discrete; k < n; ++k) {
        if (k == d.n_discrete) {
          d.c[k] = d.c[k-1] + d.p[k];
        } else {
          if (d.interpolation == Interpolation::histogram) {
            d.c[k] = d.c[k-1] + d.p[k-1]*(d.e_out[k] - d.e_out[k-1]);
          } else if (d.interpolation == Interpolation::lin_lin) {
            d.c[k] = d.c[k-1] + 0.5*(d.p[k-1] + d.p[k]) *
                  (d.e_out[k] - d.e_out[k-1]);
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
    r = (E - energy_[i]) / (energy_[i+1] - energy_[i]);
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

  // Discrete portion
  for (int j = 0; j < n_discrete; ++j) {
    k = j;
    c_k = distribution_[l].c[k];
    if (r1 < c_k) {
      end = j;
      break;
    }
  }

  // Continuous portion
  double c_k1;
  for (int j = n_discrete; j < end; ++j) {
    k = j;
    c_k1 = distribution_[l].c[k+1];
    if (r1 < c_k1) break;
    k = j + 1;
    c_k = c_k1;
  }

  double E_l_k = distribution_[l].e_out[k];
  double p_l_k = distribution_[l].p[k];
  double E_out = E_l_k;
  if (distribution_[l].interpolation == Interpolation::histogram) {
    // Histogram interpolation
    if (p_l_k > 0.0 && k >= n_discrete) {
      E_out = E_l_k + (r1 - c_k)/p_l_k;
    }

  } else if (distribution_[l].interpolation == Interpolation::lin_lin) {
    // Linear-linear interpolation
    double E_l_k1 = distribution_[l].e_out[k+1];
    double p_l_k1 = distribution_[l].p[k+1];

    if (E_l_k != E_l_k1) {
      double frac = (p_l_k1 - p_l_k)/(E_l_k1 - E_l_k);
      if (frac == 0.0) {
        E_out = E_l_k + (r1 - c_k)/p_l_k;
      } else {
        E_out = E_l_k + (std::sqrt(std::max(0.0, p_l_k*p_l_k +
                        2.0*frac*(r1 - c_k))) - p_l_k)/frac;
      }
    }
  } else {
    throw std::runtime_error{"Unexpected interpolation for continuous energy "
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
      return E_1 + (E_out - E_i_1)*(E_K - E_1)/(E_i_K - E_i_1);
    } else {
      return E_1 + (E_out - E_i1_1)*(E_K - E_1)/(E_i1_K - E_i1_1);
    }
  } else {
    return E_out;
  }
}

void ContinuousTabular::serialize(DataBuffer& buffer) const
{
  buffer.add(static_cast<int>(EnergyDistType::CONTINUOUS_TABULAR)); // 4
  buffer.add(n_region_);                                            // 4
  buffer.add(breakpoints_);                                         // 4*n_region_
  std::vector<int> interp;
  for (const auto& v : interpolation_) {
    interp.push_back(static_cast<int>(v));
  }
  buffer.add(interp);                            // 4*n_region_
  buffer.add(energy_.size());                    // 8
  buffer.add(energy_);                           // 8*energy size

  // Create locators
  std::vector<int> locators;
  int offset = 4 + 4 + (4 + 4)*n_region_ + 8 + aligned((8 + 4)*energy_.size(), 8);
  for (const auto& dist : distribution_) {
    locators.push_back(offset);
    size_t n_eout = dist.e_out.size();
    offset += 4 + 4 + 8 + 8 * 3*n_eout;
  }
  buffer.add(locators);
  buffer.align(8);

  // Write distributions
  for (const auto& dist : distribution_) {
    buffer.add(static_cast<int>(dist.interpolation));
    buffer.add(dist.n_discrete);
    buffer.add(dist.e_out.size());
    buffer.add(dist.e_out);
    buffer.add(dist.p);
    buffer.add(dist.c);
  }
}

CTTableFlat::CTTableFlat(const uint8_t* data) : data_(data)
{
  n_eout_ = *reinterpret_cast<const size_t*>(data_ + 8);
}

Interpolation CTTableFlat::interpolation() const
{
  return static_cast<Interpolation>(*reinterpret_cast<const int*>(data_));
}

int CTTableFlat::n_discrete() const
{
  return *reinterpret_cast<const int*>(data_ + 4);
}

gsl::span<const double> CTTableFlat::e_out() const
{
  auto start = reinterpret_cast<const double*>(data_ + 4 + 4 + 8);
  return {start, n_eout_};
}

gsl::span<const double> CTTableFlat::p() const
{
  auto start = reinterpret_cast<const double*>(data_ + 4 + 4 + 8 + 8*n_eout_);
  return {start, n_eout_};
}

gsl::span<const double> CTTableFlat::c() const
{
  auto start = reinterpret_cast<const double*>(data_ + 4 + 4 + 8 + (8 + 8)*n_eout_);
  return {start, n_eout_};
}

ContinuousTabularFlat::ContinuousTabularFlat(const uint8_t* data) : data_(data)
{
  n_region_ = *reinterpret_cast<const int*>(data_ + 4);
  n_energy_ = *reinterpret_cast<const size_t*>(data_ + 4 + 4 + (4 + 4)*n_region_);
}

double ContinuousTabularFlat::sample(double E, uint64_t* seed) const
{
  // Read number of interpolation regions and incoming energies
  bool histogram_interp;
  if (n_region_ == 1) {
    histogram_interp = (interpolation(0) == Interpolation::histogram);
  } else {
    histogram_interp = false;
  }

  // Find energy bin and calculate interpolation factor -- if the energy is
  // outside the range of the tabulated energies, choose the first or last bins
  auto energy_ = this->energy();
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
    r = (E - energy_[i]) / (energy_[i+1] - energy_[i]);
  }

  // Sample between the ith and [i+1]th bin
  int l;
  if (histogram_interp) {
    l = i;
  } else {
    l = r > prn(seed) ? i + 1 : i;
  }

  // Interpolation for energy E1 and EK
  auto dist_i = this->distribution(i);
  int n_discrete = dist_i.n_discrete();
  auto e_out_i = dist_i.e_out();
  int n_energy_out = e_out_i.size();
  double E_i_1 = n_energy_out > n_discrete ? e_out_i[n_discrete] : 0.0;
  double E_i_K = n_energy_out > n_discrete ? e_out_i[n_energy_out - 1] : 0.0;

  auto dist_i1 = this->distribution(i + 1);
  n_discrete = dist_i1.n_discrete();
  auto e_out_i1 = dist_i1.e_out();
  n_energy_out = e_out_i1.size();
  double E_i1_1 = n_energy_out > n_discrete ? e_out_i1[n_discrete] : 0.0;
  double E_i1_K = n_energy_out > n_discrete ? e_out_i1[n_energy_out - 1] : 0.0;

  double E_1 = E_i_1 + r * (E_i1_1 - E_i_1);
  double E_K = E_i_K + r * (E_i1_K - E_i_K);

  // Determine outgoing energy bin
  const auto& dist_l = l == i ? dist_i : dist_i1;
  auto e_out = dist_l.e_out();
  auto pdf = dist_l.p();
  auto cdf = dist_l.c();
  n_energy_out = e_out.size();
  n_discrete = dist_l.n_discrete();
  double r1 = prn(seed);
  double c_k = cdf[0];
  int k = 0;
  int end = n_energy_out - 2;

  // Discrete portion
  for (int j = 0; j < n_discrete; ++j) {
    k = j;
    c_k = cdf[k];
    if (r1 < c_k) {
      end = j;
      break;
    }
  }

  // Continuous portion
  double c_k1;
  for (int j = n_discrete; j < end; ++j) {
    k = j;
    c_k1 = cdf[k+1];
    if (r1 < c_k1) break;
    k = j + 1;
    c_k = c_k1;
  }

  double E_l_k = e_out[k];
  double p_l_k = pdf[k];
  double E_out = E_l_k;
  if (dist_l.interpolation() == Interpolation::histogram) {
    // Histogram interpolation
    if (p_l_k > 0.0 && k >= n_discrete) {
      E_out = E_l_k + (r1 - c_k)/p_l_k;
    }

  } else if (dist_l.interpolation() == Interpolation::lin_lin) {
    // Linear-linear interpolation
    double E_l_k1 = e_out[k+1];
    double p_l_k1 = pdf[k+1];

    if (E_l_k != E_l_k1) {
      double frac = (p_l_k1 - p_l_k)/(E_l_k1 - E_l_k);
      if (frac == 0.0) {
        E_out = E_l_k + (r1 - c_k)/p_l_k;
      } else {
        E_out = E_l_k + (std::sqrt(std::max(0.0, p_l_k*p_l_k +
                        2.0*frac*(r1 - c_k))) - p_l_k)/frac;
      }
    }
  } else {
    // throw std::runtime_error{"Unexpected interpolation for continuous energy "
    //   "distribution."};
  }

  // Now interpolate between incident energy bins i and i + 1
  if (!histogram_interp && n_energy_out > 1 && k >= n_discrete) {
    if (l == i) {
      return E_1 + (E_out - E_i_1)*(E_K - E_1)/(E_i_K - E_i_1);
    } else {
      return E_1 + (E_out - E_i1_1)*(E_K - E_1)/(E_i1_K - E_i1_1);
    }
  } else {
    return E_out;
  }
}

gsl::span<const int> ContinuousTabularFlat::breakpoints() const
{
  auto start = reinterpret_cast<const int*>(data_ + 4 + 4);
  return {start, n_region_};
}

Interpolation ContinuousTabularFlat::interpolation(gsl::index i) const
{
  auto start = reinterpret_cast<const int*>(data_ + 4 + 4 + 4*n_region_);
  return static_cast<Interpolation>(start[i]);
}

gsl::span<const double> ContinuousTabularFlat::energy() const
{
  auto start = reinterpret_cast<const double*>(data_ + 4 + 4 + (4 + 4)*n_region_ + 8);
  return {start, n_energy_};
}

CTTableFlat ContinuousTabularFlat::distribution(gsl::index i) const
{
  auto indices = reinterpret_cast<const int*>(data_ + 4 + 4 + (4 + 4)*n_region_ + 8 + 8*n_energy_);
  size_t offset = indices[i];
  return CTTableFlat(data_ + offset);
}

//==============================================================================
// MaxwellEnergy implementation
//==============================================================================

MaxwellEnergy::MaxwellEnergy(hid_t group)
{
  read_attribute(group, "u", u_);
  hid_t dset = open_dataset(group, "theta");
  theta_ = Tabulated1D{dset};
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
    if (E_out <= E - u_) return E_out;
  }
}

void MaxwellEnergy::serialize(DataBuffer& buffer) const
{
  buffer.add(static_cast<int>(EnergyDistType::MAXWELL)); // 4
  buffer.align(8);                                       // 4
  buffer.add(u_);                                        // 8
  theta_.serialize(buffer);
}

double MaxwellFlat::sample(double E, uint64_t* seed) const
{
  // Get temperature corresponding to incoming energy
  double theta = this->theta()(E);

  while (true) {
    // Sample maxwell fission spectrum
    double E_out = maxwell_spectrum(theta, seed);

    // Accept energy based on restriction energy
    if (E_out <= E - this->u()) return E_out;
  }
}

double MaxwellFlat::u() const
{
  return *reinterpret_cast<const double*>(data_ + 8);
}

Tabulated1DFlat MaxwellFlat::theta() const
{
  return Tabulated1DFlat(data_ + 16);
}


//==============================================================================
// Evaporation implementation
//==============================================================================

Evaporation::Evaporation(hid_t group)
{
  read_attribute(group, "u", u_);
  hid_t dset = open_dataset(group, "theta");
  theta_ = Tabulated1D{dset};
  close_dataset(dset);
}

double Evaporation::sample(double E, uint64_t* seed) const
{
  // Get temperature corresponding to incoming energy
  double theta = theta_(E);

  double y = (E - u_)/theta;
  double v = 1.0 - std::exp(-y);

  // Sample outgoing energy based on evaporation spectrum probability
  // density function
  double x;
  while (true) {
    x = -std::log((1.0 - v*prn(seed))*(1.0 - v*prn(seed)));
    if (x <= y) break;
  }

  return x * theta;
}

void Evaporation::serialize(DataBuffer& buffer) const
{
  buffer.add(static_cast<int>(EnergyDistType::EVAPORATION));  // 4
  buffer.align(8);                                            // 4
  buffer.add(u_);                                             // 8
  theta_.serialize(buffer);
}

double EvaporationFlat::sample(double E, uint64_t* seed) const
{
  // Get temperature corresponding to incoming energy
  double theta = this->theta()(E);

  double y = (E - this->u())/theta;
  double v = 1.0 - std::exp(-y);

  // Sample outgoing energy based on evaporation spectrum probability
  // density function
  double x;
  while (true) {
    x = -std::log((1.0 - v*prn(seed))*(1.0 - v*prn(seed)));
    if (x <= y) break;
  }

  return x * theta;
}

double EvaporationFlat::u() const
{
  return *reinterpret_cast<const double*>(data_ + 8);
}

Tabulated1DFlat EvaporationFlat::theta() const
{
  return Tabulated1DFlat(data_ + 16);
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
  a_ = Tabulated1D{dset};
  close_dataset(dset);
  dset = open_dataset(group, "b");
  b_ = Tabulated1D{dset};
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
    if (E_out <= E - u_) return E_out;
  }
}

void WattEnergy::serialize(DataBuffer& buffer) const
{
  buffer.add(static_cast<int>(EnergyDistType::WATT)); // 4
  buffer.align(8);                                    // 4
  buffer.add(u_);                                     // 8
  size_t n = buffer_nbytes(a_);
  buffer.add(24 + n); // offset for b                    8
  a_.serialize(buffer);
  b_.serialize(buffer);
}

double WattFlat::sample(double E, uint64_t* seed) const
{
  // Determine Watt parameters at incident energy
  double a = this->a()(E);
  double b = this->b()(E);

  double u = this->u();
  while (true) {
    // Sample energy-dependent Watt fission spectrum
    double E_out = watt_spectrum(a, b, seed);

    // Accept energy based on restriction energy
    if (E_out <= E - u) return E_out;
  }
}

double WattFlat::u() const
{
  return *reinterpret_cast<const double*>(data_ + 8);
}

Tabulated1DFlat WattFlat::a() const
{
  return Tabulated1DFlat(data_ + 24);
}

Tabulated1DFlat WattFlat::b() const
{
  auto offset = *reinterpret_cast<const size_t*>(data_ + 16);
  return Tabulated1DFlat(data_ + offset);
}

}
