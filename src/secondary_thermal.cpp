#include "openmc/secondary_thermal.h"

#include "openmc/hdf5_interface.h"
#include "openmc/random_lcg.h"
#include "openmc/search.h"
#include "openmc/vector.h"

#include "xtensor/xview.hpp"

#include <cassert>
#include <cmath> // for log, exp

namespace openmc {

// Helper function to get index on incident energy grid
void get_energy_index(
  const vector<double>& energies, double E, int& i, double& f)
{
  // Get index and interpolation factor for elastic grid
  i = 0;
  f = 0.0;
  if (E >= energies.front()) {
    i = lower_bound_index(energies.begin(), energies.end(), E);
    if (i + 1 < energies.size())
      f = (E - energies[i]) / (energies[i + 1] - energies[i]);
  }
}

double get_pdf_discrete(
  const vector<double>& mu, const vector<double>& w, double mu_0)
{
  // Make sure mu is in range [-1,1]
  if (std::abs(mu_0) > 1.0)
    mu_0 = std::copysign(1.0, mu_0);
  double a0;
  double a1;
  double b0;
  double b1;
  int32_t ai = -1;
  int32_t bi = -1;
  if (mu_0 > mu[0]) {
    ai = lower_bound_index(mu.begin(), mu.end(), mu_0);
    a0 = mu[ai];
    a1 = (ai > 1) ? mu[ai - 1] : -1.0;
  } else {
    a0 = -1.0;
    a1 = -1.0;
  }
  if (mu_0 < mu[mu.size() - 1]) {
    bi = upper_bound_index(mu.begin(), mu.end(), mu_0);
    b0 = mu[bi];
    b1 = (bi < mu.size() - 1) ? mu[bi + 1] : 1.0;
  } else {
    b0 = 1.0;
    b1 = 1.0;
  }

  //  Calculate Delta_a and Delta_b
  double delta_a = 0.5 * std::min(b0 - a0, a0 - a1);
  double delta_b = 0.5 * std::min(b1 - b0, b0 - a0);

  if (mu_0 < a0 + delta_a)
    return w[ai] / (2.0 * delta_a);
  else if (mu_0 + delta_b < b0)
    return w[bi] / (2.0 * delta_b);
  else
    return 0.0;
}

double get_pdf_discrete(const vector<double>& mu, double mu_0)
{
  vector<double> w(mu.size(), 1.0 / mu.size());
  return get_pdf_discrete(mu, w, mu_0);
}

//==============================================================================
// CoherentElasticAE implementation
//==============================================================================

CoherentElasticAE::CoherentElasticAE(const CoherentElasticXS& xs) : xs_ {xs} {}

void CoherentElasticAE::sample(
  double E_in, double& E_out, double& mu, uint64_t* seed) const
{
  // Energy doesn't change in elastic scattering (ENDF-102, Eq. 7-1)
  E_out = E_in;

  const auto& energies {xs_.bragg_edges()};

  assert(E_in >= energies.front());

  const int i = lower_bound_index(energies.begin(), energies.end(), E_in);

  // Sample a Bragg edge between 1 and i
  // E[0] < E_in < E[i+1] -> can scatter in bragg edges 0..i
  const auto& factors = xs_.factors();
  const double prob = prn(seed) * factors[i];

  const int k = std::lower_bound(factors.begin(), factors.begin() + i, prob) -
                factors.begin();

  // Characteristic scattering cosine for this Bragg edge (ENDF-102, Eq. 7-2)
  mu = 1.0 - 2.0 * energies[k] / E_in;
}

double CoherentElasticAE::get_pdf(
  double E_in, double& E_out, double& mu, uint64_t* seed) const
{
  // Energy doesn't change in elastic scattering (ENDF-102, Eq. 7-1)

  double pdf;
  E_out = E_in;
  const auto& energies {xs_.bragg_edges()};
  const auto& factors = xs_.factors();

  if (E_in < energies.front() || E_in > energies.back()) {
    return 0;
  }

  const int i = lower_bound_index(energies.begin(), energies.end(), E_in);
  vector<double> energies_cut(energies.begin(), energies.begin() + i + 1);
  vector<double> factors_cut(factors.begin(), factors.begin() + i + 1);

  vector<double> mu_vector_rev;
  std::transform(energies_cut.begin(), energies_cut.end(),
    std::back_inserter(mu_vector_rev),
    [E_in](double Ei) { return 1 - 2 * Ei / E_in; });
  vector<double> mu_vector(mu_vector_rev.rbegin(), mu_vector_rev.rend());

  auto f = xt::adapt(factors_cut, {
                                    factors_cut.size(),
                                  });
  auto weights = xt::diff(f);
  weights /= xt::sum(weights);
  vector<double> w(weights.begin(), weights.end());
  return get_pdf_discrete(mu_vector, w, mu);
}

//==============================================================================
// IncoherentElasticAE implementation
//==============================================================================

IncoherentElasticAE::IncoherentElasticAE(hid_t group)
{
  read_dataset(group, "debye_waller", debye_waller_);
}

void IncoherentElasticAE::sample(
  double E_in, double& E_out, double& mu, uint64_t* seed) const
{
  // Sample angle by inverting the distribution in ENDF-102, Eq. 7.4
  double c = 2 * E_in * debye_waller_;
  mu = std::log(1.0 + prn(seed) * (std::exp(2.0 * c) - 1)) / c - 1.0;

  // Energy doesn't change in elastic scattering (ENDF-102, Eq. 7.4)
  E_out = E_in;
}
double IncoherentElasticAE::get_pdf(
  double E_in, double& E_out, double& mu, uint64_t* seed) const
{
  // Sample angle by inverting the distribution in ENDF-102, Eq. 7.4
  double c = 2 * E_in * debye_waller_;
  E_out = E_in;

  double A = c / (1 - std::exp(-2.0 * c)); // normalization factor
  double pdf = A * std::exp(-c * (1 - mu));
  return pdf;

  // Energy doesn't change in elastic scattering (ENDF-102, Eq. 7.4)
}

//==============================================================================
// IncoherentElasticAEDiscrete implementation
//==============================================================================

IncoherentElasticAEDiscrete::IncoherentElasticAEDiscrete(
  hid_t group, const vector<double>& energy)
  : energy_ {energy}
{
  read_dataset(group, "mu_out", mu_out_);
}

void IncoherentElasticAEDiscrete::sample(
  double E_in, double& E_out, double& mu, uint64_t* seed) const
{
  // Get index and interpolation factor for elastic grid
  int i;
  double f;
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
  mu = mu_out_(i, k) + f * (mu_out_(i + 1, k) - mu_out_(i, k));

  // Inteprolate (k-1)th mu value between distributions at energies i and i+1.
  // When k==0, pick a value that will smear the cosine out to a minimum of -1.
  double mu_left = (k == 0) ? -1.0 - (mu + 1.0)
                            : mu_out_(i, k - 1) +
                                f * (mu_out_(i + 1, k - 1) - mu_out_(i, k - 1));

  // Inteprolate (k+1)th mu value between distributions at energies i and i+1.
  // When k is the last discrete value, pick a value that will smear the cosine
  // out to a maximum of 1.
  double mu_right =
    (k == n_mu - 1)
      ? 1.0 + (1.0 - mu)
      : mu_out_(i, k + 1) + f * (mu_out_(i + 1, k + 1) - mu_out_(i, k + 1));

  // Smear cosine
  mu += std::min(mu - mu_left, mu_right - mu) * (prn(seed) - 0.5);

  // Energy doesn't change in elastic scattering
  E_out = E_in;
}

double IncoherentElasticAEDiscrete::get_pdf(
  double E_in, double& E_out, double& mu, uint64_t* seed) const
{
  // Get index and interpolation factor for elastic grid
  int i;
  double f;
  get_energy_index(energy_, E_in, i, f);
  // Energy doesn't change in elastic scattering
  E_out = E_in;
  int n_mu = mu_out_.shape()[1];

  std::vector<double> mu_vector;

  for (int k = 0; k < n_mu; ++k) {
    double mu_k = mu_out_(i, k) + f * (mu_out_(i + 1, k) - mu_out_(i, k));
    mu_vector.push_back(mu_k);
  }

  return get_pdf_discrete(mu_vector, mu);
}

//==============================================================================
// IncoherentInelasticAEDiscrete implementation
//==============================================================================

IncoherentInelasticAEDiscrete::IncoherentInelasticAEDiscrete(
  hid_t group, const vector<double>& energy)
  : energy_ {energy}
{
  read_dataset(group, "energy_out", energy_out_);
  read_dataset(group, "mu_out", mu_out_);
  read_dataset(group, "skewed", skewed_);
}

void IncoherentInelasticAEDiscrete::sample(
  double E_in, double& E_out, double& mu, uint64_t* seed) const
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
    j = prn(seed) * n;
  } else {
    // Distribution skewed away from edge points
    double r = prn(seed) * (n - 3);
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
  double E_ij = energy_out_(i, j);
  double E_i1j = energy_out_(i + 1, j);

  // Outgoing energy
  E_out = (1 - f) * E_ij + f * E_i1j;

  // Sample outgoing cosine bin
  int m = mu_out_.shape()[2];
  int k = prn(seed) * m;

  // Determine outgoing cosine corresponding to E_in[i] and E_in[i+1]
  double mu_ijk = mu_out_(i, j, k);
  double mu_i1jk = mu_out_(i + 1, j, k);

  // Cosine of angle between incoming and outgoing neutron
  mu = (1 - f) * mu_ijk + f * mu_i1jk;
}

double IncoherentInelasticAEDiscrete::get_pdf(
  double E_in, double& E_out, double& mu, uint64_t* seed, int l) const
{
  // Get index and interpolation factor for inelastic grid
  int i;
  double f;
  get_energy_index(energy_, E_in, i, f);
  int j;
  int n = energy_out_.shape()[1];
  if (!skewed_) {
    // All bins equally likely
    j = prn(seed) * n;
  } else {
    // Distribution skewed away from edge points
    double r = prn(seed) * (n - 3);
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
  if (l != -1) {
    j = l;
  } // take j as param
  // Determine outgoing energy corresponding to E_in[i] and E_in[i+1]
  double E_ij = energy_out_(i, j);
  double E_i1j = energy_out_(i + 1, j);

  // Outgoing energy
  E_out = (1 - f) * E_ij + f * E_i1j;
  int m = mu_out_.shape()[2];
  std::vector<double> mu_vector;

  for (int k = 0; k < m; ++k) {
    double mu_ijk = mu_out_(i, j, k);
    double mu_i1jk = mu_out_(i + 1, j, k);
    double mu_k = (1 - f) * mu_ijk + f * mu_i1jk;
    mu_vector.push_back(mu_k);
  }

  return get_pdf_discrete(mu_vector, mu);
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

void IncoherentInelasticAE::sample(
  double E_in, double& E_out, double& mu, uint64_t* seed) const
{
  // Get index and interpolation factor for inelastic grid
  int i;
  double f;
  get_energy_index(energy_, E_in, i, f);

  // Pick closer energy based on interpolation factor
  int l = f > 0.5 ? i + 1 : i;

  // Determine outgoing energy bin
  // (First reset n_energy_out to the right value)
  auto n = distribution_[l].n_e_out;
  double r1 = prn(seed);
  double c_j = distribution_[l].e_out_cdf[0];
  double c_j1;
  std::size_t j;
  for (j = 0; j < n - 1; ++j) {
    c_j1 = distribution_[l].e_out_cdf[j + 1];
    if (r1 < c_j1)
      break;
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
    E_out = E_l_j +
            (std::sqrt(std::max(0.0, p_l_j * p_l_j + 2.0 * frac * (r1 - c_j))) -
              p_l_j) /
              frac;
  }

  // Adjustment of outgoing energy
  double E_l = energy_[l];
  if (E_out < 0.5 * E_l) {
    E_out *= 2.0 * E_in / E_l - 1.0;
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
  f = (r1 - c_j) / (c_j1 - c_j);

  // Interpolate kth mu value between distributions at energies j and j+1
  mu = mu_l(j, k) + f * (mu_l(j + 1, k) - mu_l(j, k));

  // Inteprolate (k-1)th mu value between distributions at energies j and j+1.
  // When k==0, pick a value that will smear the cosine out to a minimum of -1.
  double mu_left =
    (k == 0)
      ? mu_left = -1.0 - (mu + 1.0)
      : mu_left = mu_l(j, k - 1) + f * (mu_l(j + 1, k - 1) - mu_l(j, k - 1));

  // Inteprolate (k+1)th mu value between distributions at energies j and j+1.
  // When k is the last discrete value, pick a value that will smear the cosine
  // out to a maximum of 1.
  double mu_right =
    (k == n_mu - 1)
      ? mu_right = 1.0 + (1.0 - mu)
      : mu_right = mu_l(j, k + 1) + f * (mu_l(j + 1, k + 1) - mu_l(j, k + 1));

  // Smear cosine
  mu += std::min(mu - mu_left, mu_right - mu) * (prn(seed) - 0.5);
}

double IncoherentInelasticAE::get_pdf(
  double E_in, double& E_out, double& mu, uint64_t* seed) const
{
  // Get index and interpolation factor for inelastic grid
  int i;
  double f;
  get_energy_index(energy_, E_in, i, f);

  // Pick closer energy based on interpolation factor
  int l = f > 0.5 ? i + 1 : i;

  // Determine outgoing energy bin
  // (First reset n_energy_out to the right value)
  auto n = distribution_[l].n_e_out;
  double r1 = prn(seed);
  double c_j = distribution_[l].e_out_cdf[0];
  double c_j1;
  std::size_t j;
  for (j = 0; j < n - 1; ++j) {
    c_j1 = distribution_[l].e_out_cdf[j + 1];
    if (r1 < c_j1)
      break;
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
    E_out = E_l_j +
            (std::sqrt(std::max(0.0, p_l_j * p_l_j + 2.0 * frac * (r1 - c_j))) -
              p_l_j) /
              frac;
  }

  // Adjustment of outgoing energy
  double E_l = energy_[l];
  if (E_out < 0.5 * E_l) {
    E_out *= 2.0 * E_in / E_l - 1.0;
  } else {
    E_out += E_in - E_l;
  }

  // Sample outgoing cosine bin
  int n_mu = distribution_[l].mu.shape()[1];
  const auto& mu_l = distribution_[l].mu;
  f = (r1 - c_j) / (c_j1 - c_j);

  std::vector<double> mu_vector;

  for (int k = 0; k < n_mu; ++k) {
    double mu_k = mu_l(j, k) + f * (mu_l(j + 1, k) - mu_l(j, k));
    mu_vector.push_back(mu_k);
  }

  return get_pdf_discrete(mu_vector, mu);
}

//==============================================================================
// MixedElasticAE implementation
//==============================================================================

MixedElasticAE::MixedElasticAE(
  hid_t group, const CoherentElasticXS& coh_xs, const Function1D& incoh_xs)
  : coherent_dist_(coh_xs), coherent_xs_(coh_xs), incoherent_xs_(incoh_xs)
{
  // Read incoherent elastic distribution
  hid_t incoherent_group = open_group(group, "incoherent");
  std::string temp;
  read_attribute(incoherent_group, "type", temp);
  if (temp == "incoherent_elastic") {
    incoherent_dist_ = make_unique<IncoherentElasticAE>(incoherent_group);
  } else if (temp == "incoherent_elastic_discrete") {
    auto xs = dynamic_cast<const Tabulated1D*>(&incoh_xs);
    incoherent_dist_ =
      make_unique<IncoherentElasticAEDiscrete>(incoherent_group, xs->x());
  }
  close_group(incoherent_group);
}

void MixedElasticAE::sample(
  double E_in, double& E_out, double& mu, uint64_t* seed) const
{
  // Evaluate coherent and incoherent elastic cross sections
  double xs_coh = coherent_xs_(E_in);
  double xs_incoh = incoherent_xs_(E_in);

  if (prn(seed) * (xs_coh + xs_incoh) < xs_coh) {
    coherent_dist_.sample(E_in, E_out, mu, seed);
  } else {
    incoherent_dist_->sample(E_in, E_out, mu, seed);
  }
}

} // namespace openmc
