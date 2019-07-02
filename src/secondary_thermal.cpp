#include "openmc/secondary_thermal.h"

#include "openmc/hdf5_interface.h"
#include "openmc/random_lcg.h"
#include "openmc/search.h"

#include "xtensor/xview.hpp"

#include <cmath> // for log, exp

namespace openmc {

// Helper function to get index on incident energy grid
void
get_energy_index(const std::vector<double>& energies, double E, int& i, double& f)
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
CoherentElasticAE::sample(double E_in, double& E_out, double& mu) const
{
  // Get index and interpolation factor for elastic grid
  int i;
  double f;
  const auto& energies {xs_.bragg_edges()};
  get_energy_index(energies, E_in, i, f);

  // Sample a Bragg edge between 1 and i
  const auto& factors = xs_.factors();
  double prob = prn() * factors[i+1];
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
  read_attribute(group, "debye_waller", debye_waller_);
}

void
IncoherentElasticAE::sample(double E_in, double& E_out, double& mu) const
{
  // Sample angle by inverting the distribution in ENDF-102, Eq. 7.4
  double c = 2 * E_in * debye_waller_;
  mu = std::log(1.0 + prn()*(std::exp(2.0*c) - 1))/c - 1.0;

  // Energy doesn't change in elastic scattering (ENDF-102, Eq. 7.4)
  E_out = E_in;
}

//==============================================================================
// IncoherentElasticAEDiscrete implementation
//==============================================================================

IncoherentElasticAEDiscrete::IncoherentElasticAEDiscrete(hid_t group,
  const std::vector<double>& energy)
  : energy_{energy}
{
  read_dataset(group, "mu_out", mu_out_);
}

void
IncoherentElasticAEDiscrete::sample(double E_in, double& E_out, double& mu) const
{
  // Get index and interpolation factor for elastic grid
  int i;
  double f;
  get_energy_index(energy_, E_in, i, f);

  // Interpolate between two discrete cosines corresponding to neighboring
  // incoming energies.

  // Sample outgoing cosine bin
  int k = prn() * mu_out_.shape()[1];

  // Determine outgoing cosine corresponding to E_in[i] and E_in[i+1]
  double mu_ik  = mu_out_(i, k);
  double mu_i1k = mu_out_(i+1, k);

  // Cosine of angle between incoming and outgoing neutron
  mu = (1 - f)*mu_ik + f*mu_i1k;

  // Energy doesn't change in elastic scattering
  E_out = E_in;
}

//==============================================================================
// IncoherentInelasticAEDiscrete implementation
//==============================================================================

IncoherentInelasticAEDiscrete::IncoherentInelasticAEDiscrete(hid_t group,
  const std::vector<double>& energy)
  : energy_{energy}
{
  read_dataset(group, "energy_out", energy_out_);
  read_dataset(group, "mu_out", mu_out_);
  read_dataset(group, "skewed", skewed_);
}

void
IncoherentInelasticAEDiscrete::sample(double E_in, double& E_out, double& mu) const
{
  // Get index and interpolation factor for inelastic grid
  int i;
  double f;
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
    j = prn() * n;
  } else {
    // Distribution skewed away from edge points
    double r = prn() * (n - 3);
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
  double E_ij  = energy_out_(i, j);
  double E_i1j = energy_out_(i+1, j);

  // Outgoing energy
  E_out = (1 - f)*E_ij + f*E_i1j;

  // Sample outgoing cosine bin
  int m = mu_out_.shape()[2];
  int k = prn() * m;

  // Determine outgoing cosine corresponding to E_in[i] and E_in[i+1]
  double mu_ijk  = mu_out_(i, j, k);
  double mu_i1jk = mu_out_(i+1, j, k);

  // Cosine of angle between incoming and outgoing neutron
  mu = (1 - f)*mu_ijk + f*mu_i1jk;
}

//==============================================================================
// IncoherentInelasticAE implementation
//==============================================================================

IncoherentInelasticAE::IncoherentInelasticAE(hid_t group)
{
  // Read correlated angle-energy distribution
  CorrelatedAngleEnergy dist {group};

  // Copy incident energies
  energy_ = dist.energy();

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
          d.mu = xt::empty<double>({d.n_e_out, n_mu});
        }

        // Copy outgoing angles
        auto mu_j = xt::view(d.mu, j);
        std::copy(adist->x().begin(), adist->x().end(), mu_j.begin());
      }
    }

    distribution_.emplace_back(std::move(d));
  }

}

void
IncoherentInelasticAE::sample(double E_in, double& E_out, double& mu) const
{
  // Get index and interpolation factor for inelastic grid
  int i;
  double f;
  get_energy_index(energy_, E_in, i, f);

  // Sample between ith and [i+1]th bin
  int l = f > prn() ? i + 1 : i;

  // Determine endpoints on grid i
  auto n = distribution_[i].e_out.size();
  double E_i_1 = distribution_[i].e_out[0];
  double E_i_J = distribution_[i].e_out[n - 1];

  // Determine endpoints on grid i + 1
  n = distribution_[i + 1].e_out.size();
  double E_i1_1 = distribution_[i + 1].e_out[0];
  double E_i1_J = distribution_[i + 1].e_out[n - 1];

  double E_1 = E_i_1 + f * (E_i1_1 - E_i_1);
  double E_J = E_i_J + f * (E_i1_J - E_i_J);

  // Determine outgoing energy bin
  // (First reset n_energy_out to the right value)
  n = distribution_[l].n_e_out;
  double r1 = prn();
  double c_j = distribution_[l].e_out_cdf[0];
  double c_j1;
  std::size_t j;
  for (j = 0; j < n - 1; ++j) {
    c_j1 = distribution_[l].e_out_cdf[j + 1];
    if (r1 < c_j1) break;
    c_j = c_j1;
  }

  // check to make sure j is <= n_energy_out - 2
  j = std::min(j, n - 2);

  // Get the data to interpolate between
  double E_l_j = distribution_[l].e_out[j];
  double p_l_j = distribution_[l].e_out_pdf[j];

  // Next part assumes linear-linear interpolation in standard
  double E_l_j1 = distribution_[l].e_out[j + 1];
  double p_l_j1 = distribution_[l].e_out_pdf[j + 1];

  // Find secondary energy (variable E)
  double frac = (p_l_j1 - p_l_j) / (E_l_j1 - E_l_j);
  if (frac == 0.0) {
    E_out = E_l_j + (r1 - c_j) / p_l_j;
  } else {
    E_out = E_l_j + (std::sqrt(std::max(0.0, p_l_j*p_l_j +
          2.0*frac*(r1 - c_j))) - p_l_j) / frac;
  }

  // Now interpolate between incident energy bins i and i + 1
  if (l == i) {
    E_out = E_1 + (E_out - E_i_1) * (E_J - E_1) / (E_i_J - E_i_1);
  } else {
    E_out = E_1 + (E_out - E_i1_1) * (E_J - E_1) / (E_i1_J - E_i1_1);
  }

  // Sample outgoing cosine bin
  int n_mu = distribution_[l].mu.shape()[1];
  std::size_t k = prn() * n_mu;

  // Rather than use the sampled discrete mu directly, it is smeared over
  // a bin of width min(mu[k] - mu[k-1], mu[k+1] - mu[k]) centered on the
  // discrete mu value itself.
  const auto& mu_l = distribution_[l].mu;
  f = (r1 - c_j)/(c_j1 - c_j);

  // Determine (k-1)th mu value
  double mu_left;
  if (k == 0) {
    mu_left = -1.0;
  } else {
    mu_left = mu_l(j, k-1) + f*(mu_l(j+1, k-1) - mu_l(j, k-1));
  }

  // Determine kth mu value
  mu = mu_l(j, k) + f*(mu_l(j+1, k) - mu_l(j, k));

  // Determine (k+1)th mu value
  double mu_right;
  if (k == n_mu - 1) {
    mu_right = 1.0;
  } else {
    mu_right = mu_l(j, k+1) + f*(mu_l(j+1, k+1) - mu_l(j, k+1));
  }

  // Smear angle
  mu += std::min(mu - mu_left, mu_right - mu)*(prn() - 0.5);
}

} // namespace openmc
