#include "openmc/secondary_thermal.h"

#include "openmc/hdf5_interface.h"
#include "openmc/random_lcg.h"
#include "openmc/search.h"

#include <gsl/gsl-lite.hpp>

#include "xtensor/xview.hpp"

#include <cmath> // for log, exp

namespace openmc {

// Helper function to get index on incident energy grid
void
get_energy_index(gsl::span<const double> energies, double E, int& i, double& f)
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

//==============================================================================
// CoherentElasticAE implementation
//==============================================================================

CoherentElasticAE::CoherentElasticAE(const CoherentElasticXS& xs)
  : xs_{xs}
{ }

void
CoherentElasticAE::sample(double E_in, double& E_out, double& mu,
  uint64_t* seed) const
{
  // Energy doesn't change in elastic scattering (ENDF-102, Eq. 7-1)
  E_out = E_in;

  const auto& energies {xs_.bragg_edges()};

  Expects(E_in >= energies.front());

  const int i = lower_bound_index(energies.begin(), energies.end(), E_in);

  // Sample a Bragg edge between 1 and i
  // E[0] < E_in < E[i+1] -> can scatter in bragg edges 0..i
  const auto& factors = xs_.factors();
  const double prob = prn(seed) * factors[i];

  const int k = std::lower_bound(factors.begin(), factors.begin() + i, prob) -
                 factors.begin();

  // Characteristic scattering cosine for this Bragg edge (ENDF-102, Eq. 7-2)
  mu = 1.0 - 2.0*energies[k] / E_in;
}

void CoherentElasticAE::serialize(DataBuffer& buffer) const
{
  buffer.add(static_cast<int>(AngleEnergyType::COHERENT_ELASTIC)); // 4
  buffer.align(8);                                                 // 4
  xs_.serialize(buffer);
}

CoherentElasticXSFlat CoherentElasticAEFlat::xs() const
{
  return CoherentElasticXSFlat(data_ + 8);
}

void
CoherentElasticAEFlat::sample(double E_in, double& E_out, double& mu, uint64_t* seed) const
{
  // Get index and interpolation factor for elastic grid
  int i;
  double f;
  auto xs = this->xs();
  const auto energies = xs.bragg_edges();
  get_energy_index(energies, E_in, i, f);

  // Sample a Bragg edge between 1 and i
  const auto factors = xs.factors();
  double prob = prn(seed) * factors[i+1];
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
IncoherentElasticAE::sample(double E_in, double& E_out, double& mu,
  uint64_t* seed) const
{
  // Sample angle by inverting the distribution in ENDF-102, Eq. 7.4
  double c = 2 * E_in * debye_waller_;
  mu = std::log(1.0 + prn(seed)*(std::exp(2.0*c) - 1))/c - 1.0;

  // Energy doesn't change in elastic scattering (ENDF-102, Eq. 7.4)
  E_out = E_in;
}

void IncoherentElasticAE::serialize(DataBuffer& buffer) const
{
  buffer.add(static_cast<int>(AngleEnergyType::INCOHERENT_ELASTIC)); // 4
  buffer.align(8);                                                   // 4
  buffer.add(debye_waller_);                                         // 8
}

double IncoherentElasticAEFlat::debye_waller() const
{
  return *reinterpret_cast<const double*>(data_ + 8);
}

void
IncoherentElasticAEFlat::sample(double E_in, double& E_out, double& mu,
  uint64_t* seed) const
{
  // Sample angle by inverting the distribution in ENDF-102, Eq. 7.4
  double c = 2 * E_in * this->debye_waller();
  mu = std::log(1.0 + prn(seed)*(std::exp(2.0*c) - 1))/c - 1.0;

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
IncoherentElasticAEDiscrete::sample(double E_in, double& E_out, double& mu,
  uint64_t* seed) const
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
  mu = mu_out_(i, k) + f*(mu_out_(i+1, k) - mu_out_(i, k));

  // Inteprolate (k-1)th mu value between distributions at energies i and i+1.
  // When k==0, pick a value that will smear the cosine out to a minimum of -1.
  double mu_left = (k == 0) ?
    -1.0 - (mu + 1.0) :
    mu_out_(i, k-1) + f*(mu_out_(i+1, k-1) - mu_out_(i, k-1));

  // Inteprolate (k+1)th mu value between distributions at energies i and i+1.
  // When k is the last discrete value, pick a value that will smear the cosine
  // out to a maximum of 1.
  double mu_right = (k == n_mu - 1) ?
    1.0 + (1.0 - mu) :
    mu_out_(i, k+1) + f*(mu_out_(i+1, k+1) - mu_out_(i, k+1));

  // Smear cosine
  mu += std::min(mu - mu_left, mu_right - mu)*(prn(seed) - 0.5);

  // Energy doesn't change in elastic scattering
  E_out = E_in;
}

void IncoherentElasticAEDiscrete::serialize(DataBuffer& buffer) const
{
  buffer.add(static_cast<int>(AngleEnergyType::INCOHERENT_ELASTIC_DISCRETE)); // 4
  auto shape = mu_out_.shape();
  buffer.align(8);       // 4
  buffer.add(shape[0]);  // 8
  buffer.add(shape[1]);  // 8
  buffer.add(energy_);   // 8*...
  buffer.add(mu_out_);   // 8*...
}

IncoherentElasticAEDiscreteFlat::IncoherentElasticAEDiscreteFlat(const uint8_t* data)
  : data_(data)
{
  n_e_out_ = *reinterpret_cast<const size_t*>(data_ + 8);
  n_mu_ = *reinterpret_cast<const size_t*>(data_ + 16);
  mu_out_ = reinterpret_cast<const double*>(data_ + 24 + 8*n_e_out_);
}

gsl::span<const double> IncoherentElasticAEDiscreteFlat::energy() const
{
  auto start = reinterpret_cast<const double*>(data_ + 24);
  return {start, n_e_out_};
}

double IncoherentElasticAEDiscreteFlat::mu_out(gsl::index i, gsl::index j) const
{
  return *(mu_out_ + i*n_mu_ + j);
}

void
IncoherentElasticAEDiscreteFlat::sample(double E_in, double& E_out, double& mu,
  uint64_t* seed) const
{
  // Get index and interpolation factor for elastic grid
  int i;
  double f;
  get_energy_index(this->energy(), E_in, i, f);

  // Interpolate between two discrete cosines corresponding to neighboring
  // incoming energies.

  // Sample outgoing cosine bin
  int k = prn(seed) * n_mu_;

  // Rather than use the sampled discrete mu directly, it is smeared over
  // a bin of width 0.5*min(mu[k] - mu[k-1], mu[k+1] - mu[k]) centered on the
  // discrete mu value itself.

  // Interpolate kth mu value between distributions at energies i and i+1
  gsl::index start = 20 + 8*n_e_out_;
  mu = mu_out(i, k) + f*(mu_out(i+1, k) - mu_out(i, k));

  // Inteprolate (k-1)th mu value between distributions at energies i and i+1.
  // When k==0, pick a value that will smear the cosine out to a minimum of -1.
  double mu_left = (k == 0) ?
    -1.0 - (mu + 1.0) :
    mu_out(i, k-1) + f*(mu_out(i+1, k-1) - mu_out(i, k-1));

  // Inteprolate (k+1)th mu value between distributions at energies i and i+1.
  // When k is the last discrete value, pick a value that will smear the cosine
  // out to a maximum of 1.
  double mu_right = (k == n_mu_ - 1) ?
    1.0 + (1.0 - mu) :
    mu_out(i, k+1) + f*(mu_out(i+1, k+1) - mu_out(i, k+1));

  // Smear cosine
  mu += std::min(mu - mu_left, mu_right - mu)*(prn(seed) - 0.5);

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
IncoherentInelasticAEDiscrete::sample(double E_in, double& E_out, double& mu,
  uint64_t* seed) const
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
  double E_ij  = energy_out_(i, j);
  double E_i1j = energy_out_(i+1, j);

  // Outgoing energy
  E_out = (1 - f)*E_ij + f*E_i1j;

  // Sample outgoing cosine bin
  int m = mu_out_.shape()[2];
  int k = prn(seed) * m;

  // Determine outgoing cosine corresponding to E_in[i] and E_in[i+1]
  double mu_ijk  = mu_out_(i, j, k);
  double mu_i1jk = mu_out_(i+1, j, k);

  // Cosine of angle between incoming and outgoing neutron
  mu = (1 - f)*mu_ijk + f*mu_i1jk;
}

void IncoherentInelasticAEDiscrete::serialize(DataBuffer& buffer) const
{
  buffer.add(static_cast<int>(AngleEnergyType::INCOHERENT_INELASTIC_DISCRETE)); // 4
  auto shape = mu_out_.shape();
  buffer.add(static_cast<int>(skewed_));  // 4
  buffer.add(shape[0]);                   // 8
  buffer.add(shape[1]);                   // 8
  buffer.add(shape[2]);                   // 8
  buffer.add(energy_);                    // 8*...
  buffer.add(energy_out_);                // 8*...
  buffer.add(mu_out_);                    // 8*...
}

IncoherentInelasticAEDiscreteFlat::IncoherentInelasticAEDiscreteFlat(const uint8_t* data)
  : data_(data)
{
  n_energy_ = *reinterpret_cast<const size_t*>(data_ + 8);
  n_e_out_ = *reinterpret_cast<const size_t*>(data_ + 16);
  n_mu_ = *reinterpret_cast<const size_t*>(data_ + 24);
  energy_out_ = reinterpret_cast<const double*>(data_ + 32 + 8*n_energy_);
  mu_out_ = reinterpret_cast<const double*>(data_ + 32 + 8*n_energy_*(1 + n_e_out_));
}

bool IncoherentInelasticAEDiscreteFlat::skewed() const
{
  return *reinterpret_cast<const int*>(data_ + 4);
}

gsl::span<const double> IncoherentInelasticAEDiscreteFlat::energy() const
{
  auto start = reinterpret_cast<const double*>(data_ + 32);
  return {start, n_energy_};
}

double IncoherentInelasticAEDiscreteFlat::energy_out(gsl::index i, gsl::index j) const
{
  return *(energy_out_ + i*n_e_out_ + j);
}

double
IncoherentInelasticAEDiscreteFlat::mu_out(gsl::index i, gsl::index j, gsl::index k) const
{
  return *(mu_out_ + n_mu_*(n_e_out_*i + j) + k);
}

void
IncoherentInelasticAEDiscreteFlat::sample(double E_in, double& E_out, double& mu,
  uint64_t* seed) const
{
  // Get index and interpolation factor for inelastic grid
  int i;
  double f;
  get_energy_index(this->energy(), E_in, i, f);

  // Now that we have an incoming energy bin, we need to determine the outgoing
  // energy bin. This will depend on whether the outgoing energy distribution is
  // skewed. If it is skewed, then the first two and last two bins have lower
  // probabilities than the other bins (0.1 for the first and last bins and 0.4
  // for the second and second to last bins, relative to a normal bin
  // probability of 1). Otherwise, each bin is equally probable.

  int j;
  int n = n_e_out_;
  if (!this->skewed()) {
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
  double E_ij  = energy_out(i, j);
  double E_i1j = energy_out(i+1, j);

  // Outgoing energy
  E_out = (1 - f)*E_ij + f*E_i1j;

  // Sample outgoing cosine bin
  int m = n_mu_;
  int k = prn(seed) * m;

  // Determine outgoing cosine corresponding to E_in[i] and E_in[i+1]
  double mu_ijk  = mu_out(i, j, k);
  double mu_i1jk = mu_out(i+1, j, k);

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
IncoherentInelasticAE::sample(double E_in, double& E_out, double& mu,
  uint64_t* seed) const
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

  // Adjustment of outgoing energy
  double E_l = energy_[l];
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
  double mu_left = (k == 0) ?
    mu_left = -1.0 - (mu + 1.0) :
    mu_left = mu_l(j, k-1) + f*(mu_l(j+1, k-1) - mu_l(j, k-1));

  // Inteprolate (k+1)th mu value between distributions at energies j and j+1.
  // When k is the last discrete value, pick a value that will smear the cosine
  // out to a maximum of 1.
  double mu_right = (k == n_mu - 1) ?
    mu_right = 1.0 + (1.0 - mu) :
    mu_right = mu_l(j, k+1) + f*(mu_l(j+1, k+1) - mu_l(j, k+1));

  // Smear cosine
  mu += std::min(mu - mu_left, mu_right - mu)*(prn(seed) - 0.5);
}

void IncoherentInelasticAE::serialize(DataBuffer& buffer) const
{
  buffer.add(static_cast<int>(AngleEnergyType::INCOHERENT_INELASTIC)); // 4

  // Calculate offsets for distributions
  size_t offset = 4 + 4 + aligned((8 + 4)*energy_.size(), 8);
  std::vector<int> offsets;
  for (const auto& dist : distribution_) {
    offsets.push_back(offset);
    auto shape = dist.mu.shape();
    size_t n_e_out = shape[0];
    size_t n_mu = shape[1];
    offset += 8 + 8 + 8*n_e_out*(3 + n_mu);
  }

  // Write energy grid and offsets
  buffer.add(static_cast<int>(energy_.size()));  // 4
  buffer.add(energy_);                           // 8*energy size
  buffer.add(offsets);                           // 4*energy size
  buffer.align(8);

  // Write distributions
  for (const auto& dist : distribution_) {
    auto shape = dist.mu.shape();
    buffer.add(shape[0]);       // 8
    buffer.add(shape[1]);       // 8
    buffer.add(dist.e_out);     // 8*...
    buffer.add(dist.e_out_pdf); // 8*...
    buffer.add(dist.e_out_cdf); // 8*...
    buffer.add(dist.mu);        // 8*...
  }
}

IncoherentInelasticAEFlat::IncoherentInelasticAEFlat(const uint8_t* data)
  : data_(data)
{
  n_energy_ =  *reinterpret_cast<const size_t*>(data_ + 4);
}

gsl::span<const double> IncoherentInelasticAEFlat::energy() const
{
  auto start = reinterpret_cast<const double*>(data_ + 8);
  return {start, n_energy_};
}

DistEnergySabFlat::DistEnergySabFlat(const uint8_t* data)
  : data_(data)
{
  n_e_out_ = *reinterpret_cast<const size_t*>(data_);
  n_mu_ = *reinterpret_cast<const size_t*>(data_ + 8);
}

gsl::span<const double> DistEnergySabFlat::e_out() const
{
  auto start = reinterpret_cast<const double*>(data_ + 16);
  return {start, n_e_out_};
}

gsl::span<const double> DistEnergySabFlat::e_out_pdf() const
{
  auto start = reinterpret_cast<const double*>(data_ + 16 + 8*n_e_out_);
  return {start, n_e_out_};
}

gsl::span<const double> DistEnergySabFlat::e_out_cdf() const
{
  auto start = reinterpret_cast<const double*>(data_ + 16 + 16*n_e_out_);
  return {start, n_e_out_};
}

double DistEnergySabFlat::mu(gsl::index i, gsl::index j) const
{
  auto start = reinterpret_cast<const double*>(data_ + 16 + 24*n_e_out_);
  return *(start + i*n_mu_ + j);
}

DistEnergySabFlat IncoherentInelasticAEFlat::distribution(gsl::index i) const
{
  auto offsets = reinterpret_cast<const int*>(data_ + 8 + 8*n_energy_);
  return DistEnergySabFlat(data_ + offsets[i]);
}

void
IncoherentInelasticAEFlat::sample(double E_in, double& E_out, double& mu,
  uint64_t* seed) const
{
  // Get index and interpolation factor for inelastic grid
  int i;
  double f;
  auto energy = this->energy();
  get_energy_index(energy, E_in, i, f);

  // Pick closer energy based on interpolation factor
  int l = f > 0.5 ? i + 1 : i;

  // Determine outgoing energy bin
  // (First reset n_energy_out to the right value)
  auto dist = this->distribution(l);
  auto n = dist.n_e_out();
  double r1 = prn(seed);
  auto e_out_cdf = dist.e_out_cdf();
  double c_j = e_out_cdf[0];
  double c_j1;
  std::size_t j;
  for (j = 0; j < n - 1; ++j) {
    c_j1 = e_out_cdf[j + 1];
    if (r1 < c_j1) break;
    c_j = c_j1;
  }

  // check to make sure j is <= n_energy_out - 2
  j = std::min(j, n - 2);

  // Get the data to interpolate between
  auto e_out = dist.e_out();
  auto e_out_pdf = dist.e_out_pdf();
  double E_l_j = e_out[j];
  double p_l_j = e_out_pdf[j];

  // Next part assumes linear-linear interpolation in standard
  double E_l_j1 = e_out[j + 1];
  double p_l_j1 = e_out_pdf[j + 1];

  // Find secondary energy (variable E)
  double frac = (p_l_j1 - p_l_j) / (E_l_j1 - E_l_j);
  if (frac == 0.0) {
    E_out = E_l_j + (r1 - c_j) / p_l_j;
  } else {
    E_out = E_l_j + (std::sqrt(std::max(0.0, p_l_j*p_l_j +
          2.0*frac*(r1 - c_j))) - p_l_j) / frac;
  }

  // Adjustment of outgoing energy
  double E_l = energy[l];
  if (E_out < 0.5*E_l) {
    E_out *= 2.0*E_in/E_l - 1.0;
  } else {
    E_out += E_in - E_l;
  }

  // Sample outgoing cosine bin
  int n_mu = dist.n_mu();
  std::size_t k = prn(seed) * n_mu;

  // Rather than use the sampled discrete mu directly, it is smeared over
  // a bin of width 0.5*min(mu[k] - mu[k-1], mu[k+1] - mu[k]) centered on the
  // discrete mu value itself.
  f = (r1 - c_j)/(c_j1 - c_j);

  // Interpolate kth mu value between distributions at energies j and j+1
  mu = dist.mu(j, k) + f*(dist.mu(j+1, k) - dist.mu(j, k));

  // Inteprolate (k-1)th mu value between distributions at energies j and j+1.
  // When k==0, pick a value that will smear the cosine out to a minimum of -1.
  double mu_left = (k == 0) ?
    mu_left = -1.0 - (mu + 1.0) :
    mu_left = dist.mu(j, k-1) + f*(dist.mu(j+1, k-1) - dist.mu(j, k-1));

  // Inteprolate (k+1)th mu value between distributions at energies j and j+1.
  // When k is the last discrete value, pick a value that will smear the cosine
  // out to a maximum of 1.
  double mu_right = (k == n_mu - 1) ?
    mu_right = 1.0 + (1.0 - mu) :
    mu_right = dist.mu(j, k+1) + f*(dist.mu(j+1, k+1) - dist.mu(j, k+1));

  // Smear cosine
  mu += std::min(mu - mu_left, mu_right - mu)*(prn(seed) - 0.5);
}

} // namespace openmc
