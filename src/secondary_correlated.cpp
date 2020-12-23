#include "openmc/secondary_correlated.h"

#include <algorithm> // for copy
#include <cmath>
#include <cstddef>   // for size_t
#include <iterator>  // for back_inserter

#include <gsl/gsl>
#include "xtensor/xarray.hpp"
#include "xtensor/xview.hpp"

#include "openmc/hdf5_interface.h"
#include "openmc/distribution.h"
#include "openmc/endf.h"
#include "openmc/random_lcg.h"
#include "openmc/search.h"

namespace openmc {

//==============================================================================
//! CorrelatedAngleEnergy implementation
//==============================================================================

CorrelatedAngleEnergy::CorrelatedAngleEnergy(hid_t group)
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
  dset = open_dataset(group, "energy_out");
  std::vector<int> offsets;
  std::vector<int> interp;
  std::vector<int> n_discrete;
  read_attribute(dset, "offsets", offsets);
  read_attribute(dset, "interpolation", interp);
  read_attribute(dset, "n_discrete_lines", n_discrete);

  xt::xarray<double> eout;
  read_dataset(dset, eout);
  close_dataset(dset);

  // Read angle distributions
  xt::xarray<double> mu;
  read_dataset(group, "mu", mu);

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
    CorrTable d;
    d.interpolation = int2interp(interp[i]);
    d.n_discrete = n_discrete[i];

    // Copy data
    d.e_out = xt::view(eout, 0, xt::range(j, j+n));
    d.p = xt::view(eout, 1, xt::range(j, j+n));
    d.c = xt::view(eout, 2, xt::range(j, j+n));

    // To get answers that match ACE data, for now we still use the tabulated
    // CDF values that were passed through to the HDF5 library. At a later
    // time, we can remove the CDF values from the HDF5 library and
    // reconstruct them using the PDF
    if (false) {
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

    for (j = 0; j < n; ++j) {
      // Get interpolation scheme
      int interp_mu = std::lround(eout(3, offsets[i] + j));

      // Determine offset and size of distribution
      int offset_mu = std::lround(eout(4, offsets[i] + j));
      int m;
      if (offsets[i] + j + 1 < eout.shape()[1]) {
        m = std::lround(eout(4, offsets[i]+j+1)) - offset_mu;
      } else {
        m = mu.shape()[1] - offset_mu;
      }

      // For incoherent inelastic thermal scattering, the angle distributions
      // may be given as discrete mu values. In this case, interpolation values
      // of zero appear in the HDF5 file. Here we change it to a 1 so that
      // int2interp doesn't fail.
      if (interp_mu == 0) interp_mu = 1;

      auto interp = int2interp(interp_mu);
      auto xs = xt::view(mu, 0, xt::range(offset_mu, offset_mu + m));
      auto ps = xt::view(mu, 1, xt::range(offset_mu, offset_mu + m));
      auto cs = xt::view(mu, 2, xt::range(offset_mu, offset_mu + m));

      std::vector<double> x {xs.begin(), xs.end()};
      std::vector<double> p {ps.begin(), ps.end()};
      std::vector<double> c {cs.begin(), cs.end()};

      // To get answers that match ACE data, for now we still use the tabulated
      // CDF values that were passed through to the HDF5 library. At a later
      // time, we can remove the CDF values from the HDF5 library and
      // reconstruct them using the PDF
      Tabular* mudist = new Tabular{x.data(), p.data(), m, interp, c.data()};

      d.angle.emplace_back(mudist);
    } // outgoing energies

    distribution_.push_back(std::move(d));
  } // incoming energies
}

void CorrelatedAngleEnergy::sample(double E_in, double& E_out, double& mu,
  uint64_t* seed) const
{
  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< REMOVE THIS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  // Before the secondary distribution refactor, an isotropic polar cosine was
  // always sampled but then overwritten with the polar cosine sampled from the
  // correlated distribution. To preserve the random number stream, we keep
  // this dummy sampling here but can remove it later (will change answers)
  mu = 2.0*prn(seed) - 1.0;
  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< REMOVE THIS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  // Find energy bin and calculate interpolation factor -- if the energy is
  // outside the range of the tabulated energies, choose the first or last bins
  auto n_energy_in = energy_.size();
  int i;
  double r;
  if (E_in < energy_[0]) {
    i = 0;
    r = 0.0;
  } else if (E_in > energy_[n_energy_in - 1]) {
    i = n_energy_in - 2;
    r = 1.0;
  } else {
    i = lower_bound_index(energy_.begin(), energy_.end(), E_in);
    r = (E_in - energy_[i]) / (energy_[i+1] - energy_[i]);
  }

  // Sample between the ith and [i+1]th bin
  int l = r > prn(seed) ? i + 1 : i;

  // Interpolation for energy E1 and EK
  int n_energy_out = distribution_[i].e_out.size();
  int n_discrete = distribution_[i].n_discrete;
  double E_i_1 = distribution_[i].e_out[n_discrete];
  double E_i_K = distribution_[i].e_out[n_energy_out - 1];

  n_energy_out = distribution_[i+1].e_out.size();
  n_discrete = distribution_[i+1].n_discrete;
  double E_i1_1 = distribution_[i+1].e_out[n_discrete];
  double E_i1_K = distribution_[i+1].e_out[n_energy_out - 1];

  double E_1 = E_i_1 + r*(E_i1_1 - E_i_1);
  double E_K = E_i_K + r*(E_i1_K - E_i_K);

  // Determine outgoing energy bin
  n_energy_out = distribution_[l].e_out.size();
  n_discrete = distribution_[l].n_discrete;
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
  if (distribution_[l].interpolation == Interpolation::histogram) {
    // Histogram interpolation
    if (p_l_k > 0.0 && k >= n_discrete) {
      E_out = E_l_k + (r1 - c_k)/p_l_k;
    } else {
      E_out = E_l_k;
    }

  } else if (distribution_[l].interpolation == Interpolation::lin_lin) {
    // Linear-linear interpolation
    double E_l_k1 = distribution_[l].e_out[k+1];
    double p_l_k1 = distribution_[l].p[k+1];

    double frac = (p_l_k1 - p_l_k)/(E_l_k1 - E_l_k);
    if (frac == 0.0) {
      E_out = E_l_k + (r1 - c_k)/p_l_k;
    } else {
      E_out = E_l_k + (std::sqrt(std::max(0.0, p_l_k*p_l_k +
                        2.0*frac*(r1 - c_k))) - p_l_k)/frac;
    }

  }

  // Now interpolate between incident energy bins i and i + 1
  if (k >= n_discrete){
    if (l == i) {
      E_out = E_1 + (E_out - E_i_1)*(E_K - E_1)/(E_i_K - E_i_1);
    } else {
      E_out = E_1 + (E_out - E_i1_1)*(E_K - E_1)/(E_i1_K - E_i1_1);
    }
  }

  // Find correlated angular distribution for closest outgoing energy bin
  if (r1 - c_k < c_k1 - r1) {
    mu = distribution_[l].angle[k]->sample(seed);
  } else {
    mu = distribution_[l].angle[k + 1]->sample(seed);
  }
}

UnifiedAngleEnergy CorrelatedAngleEnergy::serialize() const
{
  // Determine size of buffer needed
  size_t n = 4 + (4 + 4)*n_region_ + 8 + (8 + 4)*energy_.size();
  int offset = n;
  std::vector<int> locators;
  for (const auto& dist : distribution_) {
    locators.push_back(n);
    size_t n_eout = dist.e_out.size();
    n += 4 + 4 + 8 + (8*3 + 4)*n_eout;

    for (const auto& adist : dist.angle) {
      n += adist->nbytes();
    };
  }
  DataBuffer buffer(n);

  // Write interpolation information
  buffer.add(n_region_);
  buffer.add(breakpoints_);
  std::vector<int> interp;
  for (auto v : interpolation_) {
    interp.push_back(static_cast<int>(v));
  }
  buffer.add(interp);

  // Write incident energies and locators
  buffer.add(energy_.size());
  buffer.add(energy_);
  buffer.add(locators);

  // Write distributions
  for (const auto& dist : distribution_) {
    buffer.add(dist.n_discrete);
    buffer.add(static_cast<int>(dist.interpolation));
    buffer.add(dist.e_out.size());
    buffer.add(dist.e_out);
    buffer.add(dist.p);
    buffer.add(dist.c);

    // Write locators for angle distributions
    std::vector<int> angle_dist_offsets;
    int offset = 4 + 4 + 8 + (8*3 + 4)*dist.e_out.size();
    for (const auto& adist : dist.angle) {
      angle_dist_offsets.push_back(offset);
      offset += adist->nbytes();
    }
    buffer.add(angle_dist_offsets);

    // Write angle distributions
    for (const auto& adist : dist.angle) {
      adist->serialize(buffer);
    }
  }
  Ensures(n == buffer.offset_);

  return {AngleEnergyType::CORRELATED, std::move(buffer)};
}

CorrTableFlat::CorrTableFlat(const uint8_t* data) : data_(data)
{
  n_eout_ = *reinterpret_cast<const size_t*>(data_ + 8);
}

int CorrTableFlat::n_discrete() const
{
  return *reinterpret_cast<const int*>(data_);
}

Interpolation CorrTableFlat::interpolation() const
{
  return static_cast<Interpolation>(*reinterpret_cast<const int*>(data_ + 4));
}

gsl::span<const double> CorrTableFlat::e_out() const
{
  auto start = reinterpret_cast<const double*>(data_ + 4 + 4 + 8);
  return {start, n_eout_};
}

gsl::span<const double> CorrTableFlat::p() const
{
  auto start = reinterpret_cast<const double*>(data_ + 4 + 4 + 8 + 8*n_eout_);
  return {start, n_eout_};
}

gsl::span<const double> CorrTableFlat::c() const
{
  auto start = reinterpret_cast<const double*>(data_ + 4 + 4 + 8 + 16*n_eout_);
  return {start, n_eout_};
}

TabularFlat CorrTableFlat::angle(gsl::index i) const
{
  auto indices = reinterpret_cast<const int*>(data_ + 4 + 4 + 8 + 24*n_eout_);
  size_t offset = indices[i];
  return TabularFlat(data_ + offset);
}

CorrelatedAngleEnergyFlat::CorrelatedAngleEnergyFlat(const uint8_t* data) : data_(data)
{
  n_region_ = *reinterpret_cast<const int*>(data_);
  n_energy_ = *reinterpret_cast<const size_t*>(data_ + 4 + (4 + 4)*n_region_);
}

void CorrelatedAngleEnergyFlat::sample(double E_in, double& E_out, double& mu, uint64_t* seed) const
{
  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< REMOVE THIS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  // Before the secondary distribution refactor, an isotropic polar cosine was
  // always sampled but then overwritten with the polar cosine sampled from the
  // correlated distribution. To preserve the random number stream, we keep
  // this dummy sampling here but can remove it later (will change answers)
  mu = 2.0*prn(seed) - 1.0;
  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< REMOVE THIS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  // Find energy bin and calculate interpolation factor -- if the energy is
  // outside the range of the tabulated energies, choose the first or last bins
  auto energy_ = this->energy();
  auto n_energy_in = energy_.size();
  int i;
  double r;
  if (E_in < energy_[0]) {
    i = 0;
    r = 0.0;
  } else if (E_in > energy_[n_energy_in - 1]) {
    i = n_energy_in - 2;
    r = 1.0;
  } else {
    i = lower_bound_index(energy_.begin(), energy_.end(), E_in);
    r = (E_in - energy_[i]) / (energy_[i+1] - energy_[i]);
  }

  // Sample between the ith and [i+1]th bin
  int l = r > prn(seed) ? i + 1 : i;

  // Interpolation for energy E1 and EK
  auto dist_i = this->distribution(i);
  int n_discrete = dist_i.n_discrete();
  auto e_out_i = dist_i.e_out();
  int n_energy_out = e_out_i.size();
  double E_i_1 = e_out_i[n_discrete];
  double E_i_K = e_out_i[n_energy_out - 1];

  auto dist_i1 = this->distribution(i + 1);
  n_discrete = dist_i1.n_discrete();
  auto e_out_i1 = dist_i1.e_out();
  n_energy_out = e_out_i1.size();
  double E_i1_1 = e_out_i1[n_discrete];
  double E_i1_K = e_out_i1[n_energy_out - 1];

  double E_1 = E_i_1 + r*(E_i1_1 - E_i_1);
  double E_K = E_i_K + r*(E_i1_K - E_i_K);

  // Determine outgoing energy bin
  const auto& dist_l = l == i ? dist_i : dist_i1;
  auto e_out = dist_l.e_out();
  auto pdf = dist_l.p();
  auto cdf = dist_l.c();
  n_discrete = dist_l.n_discrete();
  n_energy_out = e_out.size();
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
  if (dist_l.interpolation() == Interpolation::histogram) {
    // Histogram interpolation
    if (p_l_k > 0.0 && k >= n_discrete) {
      E_out = E_l_k + (r1 - c_k)/p_l_k;
    } else {
      E_out = E_l_k;
    }

  } else if (dist_l.interpolation() == Interpolation::lin_lin) {
    // Linear-linear interpolation
    double E_l_k1 = e_out[k+1];
    double p_l_k1 = pdf[k+1];

    double frac = (p_l_k1 - p_l_k)/(E_l_k1 - E_l_k);
    if (frac == 0.0) {
      E_out = E_l_k + (r1 - c_k)/p_l_k;
    } else {
      E_out = E_l_k + (std::sqrt(std::max(0.0, p_l_k*p_l_k +
                        2.0*frac*(r1 - c_k))) - p_l_k)/frac;
    }

  }

  // Now interpolate between incident energy bins i and i + 1
  if (k >= n_discrete){
    if (l == i) {
      E_out = E_1 + (E_out - E_i_1)*(E_K - E_1)/(E_i_K - E_i_1);
    } else {
      E_out = E_1 + (E_out - E_i1_1)*(E_K - E_1)/(E_i1_K - E_i1_1);
    }
  }

  // Find correlated angular distribution for closest outgoing energy bin
  if (r1 - c_k < c_k1 - r1) {
    mu = dist_l.angle(k).sample(seed);
  } else {
    mu = dist_l.angle(k + 1).sample(seed);
  }
}

gsl::span<const int> CorrelatedAngleEnergyFlat::breakpoints() const
{
  auto start = reinterpret_cast<const int*>(data_ + 4);
  return {start, n_region_};
}

Interpolation CorrelatedAngleEnergyFlat::interpolation(gsl::index i) const
{
  auto start = reinterpret_cast<const int*>(data_ + 4 + 4*n_region_);
  return static_cast<Interpolation>(start[i]);
}

gsl::span<const double> CorrelatedAngleEnergyFlat::energy() const
{
  auto start = reinterpret_cast<const double*>(data_ + 4 + (4 + 4)*n_region_ + 8);
  return {start, n_energy_};
}

CorrTableFlat CorrelatedAngleEnergyFlat::distribution(gsl::index i) const
{
  auto indices = reinterpret_cast<const int*>(data_ + 4 + (4 + 4)*n_region_ + 8 + 8*n_energy_);
  size_t offset = indices[i];
  return CorrTableFlat(data_ + offset);
}

} // namespace openmc
