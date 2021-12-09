#include "openmc/secondary_thermal.h"

#include "openmc/hdf5_interface.h"
#include "openmc/random_lcg.h"
#include "openmc/search.h"

#include "xtensor/xview.hpp"

#include <cmath> // for log, exp

namespace openmc {

// Helper function to get index on incident energy grid
HD void get_energy_index(
  const vector<xsfloat>& energies, xsfloat E, int& i, xsfloat& f)
{
  // Get index and interpolation factor for elastic grid
  i = 0;
  f = 0.0;
  if (E >= energies.front()) {
    i = lower_bound_index(energies.begin(), energies.end(), E);
    f = (E - energies[i]) / (energies[i+1] - energies[i]);
  }
}

//==============================================================================
// CoherentElasticAE implementation
//==============================================================================

CoherentElasticAE::CoherentElasticAE(const CoherentElasticXS& xs)
  : xs_{xs}
{ }

void
CoherentElasticAE::sample(xsfloat E_in, xsfloat& E_out, xsfloat& mu,
  uint64_t* seed) const
{
  // Get index and interpolation factor for elastic grid
  int i;
  xsfloat f;
  const auto& energies {xs_.bragg_edges()};
  get_energy_index(energies, E_in, i, f);

  // Sample a Bragg edge between 1 and i
  const auto& factors = xs_.factors();
  xsfloat prob = prn(seed) * factors[i+1];
  int k = 0;
  if (prob >= factors.front()) {
    k = lower_bound_index(factors.begin(), factors.begin() + (i+1), prob);
  }

  // Characteristic scattering cosine for this Bragg edge (ENDF-102, Eq. 7-2)
  mu = 1.0 - 2.0*energies[k] / E_in;

  // Energy doesn't change in elastic scattering (ENDF-102, Eq. 7-1)
  E_out = E_in;
}

//==============================================================================
// IncoherentElasticAE implementation
//==============================================================================

IncoherentElasticAE::IncoherentElasticAE(hid_t group)
{
  read_dataset(group, "debye_waller", debye_waller_);
}

void
IncoherentElasticAE::sample(xsfloat E_in, xsfloat& E_out, xsfloat& mu,
  uint64_t* seed) const
{
  // Sample angle by inverting the distribution in ENDF-102, Eq. 7.4
  xsfloat c = 2 * E_in * debye_waller_;
  mu = std::log(1.0 + prn(seed)*(std::exp(2.0*c) - 1))/c - 1.0;

  // Energy doesn't change in elastic scattering (ENDF-102, Eq. 7.4)
  E_out = E_in;
}

//==============================================================================
// IncoherentElasticAEDiscrete implementation
//==============================================================================

IncoherentElasticAEDiscrete::IncoherentElasticAEDiscrete(
  hid_t group, const vector<xsfloat>& energy)
  : energy_ {energy}
{
  read_dataset(group, "mu_out", mu_out_);
}

void
IncoherentElasticAEDiscrete::sample(xsfloat E_in, xsfloat& E_out, xsfloat& mu,
  uint64_t* seed) const
{
  // Get index and interpolation factor for elastic grid
  int i;
  xsfloat f;
  get_energy_index(energy_, E_in, i, f);

  // Interpolate between two discrete cosines corresponding to neighboring
  // incoming energies.

  // Sample outgoing cosine bin
  int n_mu = mu_out_.shape()[1];
  int k = prn(seed) * n_mu;

  // Rather than use the sampled discrete mu directly, it is smeared over
  // a bin of width 0.5*min(mu[k] - mu[k-1], mu[k+1] - mu[k]) centered on the
  // discrete mu value itself.

  // Interpolate kth mu value between distributions at energies i and i+1
  mu = mu_out_(i, k) + f*(mu_out_(i+1, k) - mu_out_(i, k));

  // Inteprolate (k-1)th mu value between distributions at energies i and i+1.
  // When k==0, pick a value that will smear the cosine out to a minimum of -1.
  xsfloat mu_left = (k == 0) ?
    -1.0 - (mu + 1.0) :
    mu_out_(i, k-1) + f*(mu_out_(i+1, k-1) - mu_out_(i, k-1));

  // Inteprolate (k+1)th mu value between distributions at energies i and i+1.
  // When k is the last discrete value, pick a value that will smear the cosine
  // out to a maximum of 1.
  xsfloat mu_right = (k == n_mu - 1) ?
    1.0 + (1.0 - mu) :
    mu_out_(i, k+1) + f*(mu_out_(i+1, k+1) - mu_out_(i, k+1));

  // Smear cosine
  mu += std::min(mu - mu_left, mu_right - mu)*(prn(seed) - 0.5);

  // Energy doesn't change in elastic scattering
  E_out = E_in;
}

//==============================================================================
// IncoherentInelasticAEDiscrete implementation
//==============================================================================

IncoherentInelasticAEDiscrete::IncoherentInelasticAEDiscrete(
  hid_t group, const vector<xsfloat>& energy)
  : energy_ {energy}
{
  read_dataset(group, "energy_out", energy_out_);
  read_dataset(group, "mu_out", mu_out_);
  read_dataset(group, "skewed", skewed_);
}

void
IncoherentInelasticAEDiscrete::sample(xsfloat E_in, xsfloat& E_out, xsfloat& mu,
  uint64_t* seed) const
{
  // Get index and interpolation factor for inelastic grid
  int i;
  xsfloat f;
  get_energy_index(energy_, E_in, i, f);

  // Now that we have an incoming energy bin, we need to determine the outgoing
  // energy bin. This will depend on whether the outgoing energy distribution is
  // skewed. If it is skewed, then the first two and last two bins have lower
  // probabilities than the other bins (0.1 for the first and last bins and 0.4
  // for the second and second to last bins, relative to a normal bin
  // probability of 1). Otherwise, each bin is equally probable.

  int j;
  int n = energy_out_.shape()[1];
  if (!skewed_) {
    // All bins equally likely
    j = prn(seed) * n;
  } else {
    // Distribution skewed away from edge points
    xsfloat r = prn(seed) * (n - 3);
    if (r > 1.0) {
      // equally likely N-4 middle bins
      j = r + 1;
    } else if (r > 0.6) {
      // second to last bin has relative probability of 0.4
      j = n - 2;
    } else if (r > 0.5) {
      // last bin has relative probability of 0.1
      j = n - 1;
    } else if (r > 0.1) {
      // second bin has relative probability of 0.4
      j = 1;
    } else {
      // first bin has relative probability of 0.1
      j = 0;
    }
  }

  // Determine outgoing energy corresponding to E_in[i] and E_in[i+1]
  xsfloat E_ij  = energy_out_(i, j);
  xsfloat E_i1j = energy_out_(i+1, j);

  // Outgoing energy
  E_out = (1 - f)*E_ij + f*E_i1j;

  // Sample outgoing cosine bin
  int m = mu_out_.shape()[2];
  int k = prn(seed) * m;

  // Determine outgoing cosine corresponding to E_in[i] and E_in[i+1]
  xsfloat mu_ijk  = mu_out_(i, j, k);
  xsfloat mu_i1jk = mu_out_(i+1, j, k);

  // Cosine of angle between incoming and outgoing neutron
  mu = (1 - f)*mu_ijk + f*mu_i1jk;
}

//==============================================================================
// IncoherentInelasticAE implementation
//==============================================================================

IncoherentInelasticAE::IncoherentInelasticAE(hid_t group)
{
  hid_t dset = open_dataset(group, "energy");
  read_dataset(dset, energy_);
  close_dataset(dset);

  dset = open_dataset(group, "energy_out");
  vector<int> offsets;
  read_attribute(dset, "offsets", offsets);
  xt::xarray<xsfloat> eout;
  read_dataset(dset, eout);
  close_dataset(dset);

  xt::xarray<xsfloat> mu;
  read_dataset(group, "mu", mu);

  distribution_.reserve(energy_.size());

  auto n_energy = energy_.size();
  auto n_energy_out = eout.shape()[1];
  for (int i = 0; i < n_energy; ++i) {
    // Determine number of outgoing energies
    int j = offsets[i];
    int n;
    if (i < n_energy - 1) {
      n = offsets[i + 1] - j;
    } else {
      n = eout.shape()[1] - j;
    }
    DistEnergySab d;

    d.n_e_out = n;
    fill_vec_from_view(d.e_out, xt::view(eout, 0, xt::range(j, j + n)));
    fill_vec_from_view(d.e_out_pdf, xt::view(eout, 1, xt::range(j, j + n)));
    fill_vec_from_view(d.e_out_cdf, xt::view(eout, 2, xt::range(j, j + n)));
    for (decltype(d.e_out_pdf)::size_type ii = 0; ii < d.e_out_pdf.size();
         ++ii) {
      d.e_out_pdf[ii] /= d.e_out_cdf[n - 1];
      d.e_out_cdf[ii] /= d.e_out_cdf[n - 1];
    }

    for (j = 0; j < n; ++j) { // Loop over e_out  values
      int offset_mu = std::lround(eout(4, offsets[i] + j));
      int m; // number of mu values
      if (offsets[i] + j + 1 < eout.shape()[1]) {
        m = std::lround(eout(4, offsets[i] + j + 1)) - offset_mu;
      } else {
        m = mu.shape()[1] - offset_mu;
      }

      if (j == 0) {
        std::array<decltype(d.mu)::size_type, 2> mu_shape = {
          static_cast<unsigned>(d.n_e_out), static_cast<unsigned>(m)};
        d.mu.resize(mu_shape);
      }

      // Now fill mu values into the tensor
      auto shape1 = mu.shape();
      auto shape2 = d.mu.shape();
      for (int i_mu = 0; i_mu < m; ++i_mu) {
        d.mu(j, i_mu) = mu(0, offset_mu + i_mu);
      }
    }
    distribution_.push_back(std::move(d));
  }
}

void
IncoherentInelasticAE::sample(xsfloat E_in, xsfloat& E_out, xsfloat& mu,
  uint64_t* seed) const
{
  // Get index and interpolation factor for inelastic grid
  int i;
  xsfloat f;
  get_energy_index(energy_, E_in, i, f);

  // Pick closer energy based on interpolation factor
  int l = f > 0.5 ? i + 1 : i;

  // Determine outgoing energy bin
  // (First reset n_energy_out to the right value)
  auto n = distribution_[l].n_e_out;
  xsfloat r1 = prn(seed);
  xsfloat c_j = distribution_[l].e_out_cdf[0];
  xsfloat c_j1;
  std::size_t j;
  for (j = 0; j < n - 1; ++j) {
    c_j1 = distribution_[l].e_out_cdf[j + 1];
    if (r1 < c_j1) break;
    c_j = c_j1;
  }

  // check to make sure j is <= n_energy_out - 2
  j = std::min(j, n - 2);

  // Get the data to interpolate between
  xsfloat E_l_j = distribution_[l].e_out[j];
  xsfloat p_l_j = distribution_[l].e_out_pdf[j];

  // Next part assumes linear-linear interpolation in standard
  xsfloat E_l_j1 = distribution_[l].e_out[j + 1];
  xsfloat p_l_j1 = distribution_[l].e_out_pdf[j + 1];

  // Find secondary energy (variable E)
  xsfloat frac = (p_l_j1 - p_l_j) / (E_l_j1 - E_l_j);
  if (frac == 0.0) {
    E_out = E_l_j + (r1 - c_j) / p_l_j;
  } else {
    E_out = E_l_j + (std::sqrt(std::max(0.0, p_l_j*p_l_j +
          2.0*frac*(r1 - c_j))) - p_l_j) / frac;
  }

  // Adjustment of outgoing energy
  xsfloat E_l = energy_[l];
  if (E_out < 0.5*E_l) {
    E_out *= 2.0*E_in/E_l - 1.0;
  } else {
    E_out += E_in - E_l;
  }

  // Sample outgoing cosine bin
  int n_mu = distribution_[l].mu.shape()[1];
  std::size_t k = prn(seed) * n_mu;

  // Rather than use the sampled discrete mu directly, it is smeared over
  // a bin of width 0.5*min(mu[k] - mu[k-1], mu[k+1] - mu[k]) centered on the
  // discrete mu value itself.
  const auto& mu_l = distribution_[l].mu;
  f = (r1 - c_j)/(c_j1 - c_j);

  // Interpolate kth mu value between distributions at energies j and j+1
  mu = mu_l(j, k) + f*(mu_l(j+1, k) - mu_l(j, k));

  // Inteprolate (k-1)th mu value between distributions at energies j and j+1.
  // When k==0, pick a value that will smear the cosine out to a minimum of -1.
  xsfloat mu_left = (k == 0) ?
    mu_left = -1.0 - (mu + 1.0) :
    mu_left = mu_l(j, k-1) + f*(mu_l(j+1, k-1) - mu_l(j, k-1));

  // Inteprolate (k+1)th mu value between distributions at energies j and j+1.
  // When k is the last discrete value, pick a value that will smear the cosine
  // out to a maximum of 1.
  xsfloat mu_right = (k == n_mu - 1) ?
    mu_right = 1.0 + (1.0 - mu) :
    mu_right = mu_l(j, k+1) + f*(mu_l(j+1, k+1) - mu_l(j, k+1));

  // Smear cosine
  mu += std::min(mu - mu_left, mu_right - mu)*(prn(seed) - 0.5);
}

} // namespace openmc
