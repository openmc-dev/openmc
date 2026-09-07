//! \file secondary_thermal.h
//! Angle-energy distributions for thermal scattering

#ifndef OPENMC_SECONDARY_THERMAL_H
#define OPENMC_SECONDARY_THERMAL_H

#include "openmc/angle_energy.h"
#include "openmc/endf.h"
#include "openmc/search.h"
#include "openmc/secondary_correlated.h"
#include "openmc/vector.h"

#include "openmc/tensor.h"
#include <hdf5.h>

namespace openmc {

//==============================================================================
//! Coherent elastic scattering angle-energy distribution
//==============================================================================

class CoherentElasticAE : public AngleEnergy {
public:
  //! Construct from a coherent elastic scattering cross section
  //
  //! \param[in] xs Coherent elastic scattering cross section
  explicit CoherentElasticAE(const CoherentElasticXS& xs);

  //! Sample distribution for an angle and energy
  //! \param[in] E_in Incoming energy in [eV]
  //! \param[out] E_out Outgoing energy in [eV]
  //! \param[out] mu Outgoing cosine with respect to current direction
  //! \param[inout] seed Pseudorandom seed pointer
  void sample(
    double E_in, double& E_out, double& mu, uint64_t* seed) const override;

  //! Sample an outgoing energy and evaluate the angular PDF
  //! \param[in] E_in Incoming energy in [eV]
  //! \param[in] mu Scattering cosine with respect to current direction
  //! \param[out] E_out Outgoing energy in [eV]
  //! \param[inout] seed Pseudorandom seed pointer
  //! \return Probability density for the scattering cosine
  double sample_energy_and_pdf(
    double E_in, double mu, double& E_out, uint64_t* seed) const override;

private:
  const CoherentElasticXS& xs_; //!< Coherent elastic scattering cross section
  tensor::Tensor<double> bragg_edges_; //!< Copy of Bragg edges for slicing
  tensor::Tensor<double>
    factors_diff_; //!< Differences over elastic scattering factors
};

//==============================================================================
//! Incoherent elastic scattering angle-energy distribution
//==============================================================================

class IncoherentElasticAE : public AngleEnergy {
public:
  //! Construct from HDF5 file
  //
  //! \param[in] group  HDF5 group
  explicit IncoherentElasticAE(hid_t group);

  //! Sample distribution for an angle and energy
  //! \param[in] E_in Incoming energy in [eV]
  //! \param[out] E_out Outgoing energy in [eV]
  //! \param[out] mu Outgoing cosine with respect to current direction
  //! \param[inout] seed Pseudorandom number seed pointer
  void sample(
    double E_in, double& E_out, double& mu, uint64_t* seed) const override;

  //! Sample an outgoing energy and evaluate the angular PDF
  //! \param[in] E_in Incoming energy in [eV]
  //! \param[in] mu Scattering cosine with respect to current direction
  //! \param[out] E_out Outgoing energy in [eV]
  //! \param[inout] seed Pseudorandom seed pointer
  //! \return Probability density for the scattering cosine
  double sample_energy_and_pdf(
    double E_in, double mu, double& E_out, uint64_t* seed) const override;

private:
  double debye_waller_;
};

//==============================================================================
//! Incoherent elastic scattering angle-energy distribution (discrete)
//==============================================================================

class IncoherentElasticAEDiscrete : public AngleEnergy {
public:
  //! Construct from HDF5 file
  //
  //! \param[in] group  HDF5 group
  //! \param[in] energy  Energies at which cosines are tabulated
  explicit IncoherentElasticAEDiscrete(
    hid_t group, const vector<double>& energy);

  //! Sample distribution for an angle and energy
  //! \param[in] E_in Incoming energy in [eV]
  //! \param[out] E_out Outgoing energy in [eV]
  //! \param[out] mu Outgoing cosine with respect to current direction
  //! \param[inout] seed Pseudorandom number seed pointer
  void sample(
    double E_in, double& E_out, double& mu, uint64_t* seed) const override;

  //! Sample an outgoing energy and evaluate the angular PDF
  //! \param[in] E_in Incoming energy in [eV]
  //! \param[in] mu Scattering cosine with respect to current direction
  //! \param[out] E_out Outgoing energy in [eV]
  //! \param[inout] seed Pseudorandom seed pointer
  //! \return Probability density for the scattering cosine
  double sample_energy_and_pdf(
    double E_in, double mu, double& E_out, uint64_t* seed) const override;

private:
  const vector<double>& energy_;  //!< Energies at which cosines are tabulated
  tensor::Tensor<double> mu_out_; //!< Cosines for each incident energy
};

//==============================================================================
//! Incoherent inelastic scattering angle-energy distribution (discrete)
//==============================================================================

class IncoherentInelasticAEDiscrete : public AngleEnergy {
public:
  //! Construct from HDF5 file
  //
  //! \param[in] group  HDF5 group
  //! \param[in] energy  Incident energies at which distributions are tabulated
  explicit IncoherentInelasticAEDiscrete(
    hid_t group, const vector<double>& energy);

  //! Sample distribution for an angle and energy
  //! \param[in] E_in Incoming energy in [eV]
  //! \param[out] E_out Outgoing energy in [eV]
  //! \param[out] mu Outgoing cosine with respect to current direction
  //! \param[inout] seed Pseudorandom number seed pointer
  void sample(
    double E_in, double& E_out, double& mu, uint64_t* seed) const override;
  //! Sample outgoing energy bin parameters
  //! \param[in] E_in Incoming energy in [eV]
  //! \param[out] E_out Outgoing energy in [eV]
  //! \param[out] j Sampled outgoing energy bin index
  //! \param[inout] seed Pseudorandom seed pointer
  void sample_params(double E_in, double& E_out, int& j, uint64_t* seed) const;

  //! Sample an outgoing energy and evaluate the angular PDF
  //! \param[in] E_in Incoming energy in [eV]
  //! \param[in] mu Scattering cosine with respect to current direction
  //! \param[out] E_out Outgoing energy in [eV]
  //! \param[inout] seed Pseudorandom seed pointer
  //! \return Probability density for the scattering cosine
  double sample_energy_and_pdf(
    double E_in, double mu, double& E_out, uint64_t* seed) const override;

private:
  const vector<double>& energy_; //!< Incident energies
  tensor::Tensor<double>
    energy_out_; //!< Outgoing energies for each incident energy
  tensor::Tensor<double>
    mu_out_;    //!< Outgoing cosines for each incident/outgoing energy
  bool skewed_; //!< Whether outgoing energy distribution is skewed
};

//==============================================================================
//! Incoherent inelastic scattering angle-energy distribution
//==============================================================================

class IncoherentInelasticAE : public AngleEnergy {
public:
  //! Construct from HDF5 file
  //
  //! \param[in] group  HDF5 group
  explicit IncoherentInelasticAE(hid_t group);

  //! Sample distribution for an angle and energy
  //! \param[in] E_in Incoming energy in [eV]
  //! \param[out] E_out Outgoing energy in [eV]
  //! \param[out] mu Outgoing cosine with respect to current direction
  //! \param[inout] seed Pseudorandom number seed pointer
  void sample(
    double E_in, double& E_out, double& mu, uint64_t* seed) const override;

  //! Sample outgoing energy bin parameters
  //! \param[in] E_in Incoming energy in [eV]
  //! \param[out] E_out Outgoing energy in [eV]
  //! \param[out] f Interpolation factor within sampled energy bin
  //! \param[out] l Index of the closer incident energy
  //! \param[out] j Sampled outgoing energy bin index
  //! \param[inout] seed Pseudorandom seed pointer
  void sample_params(double E_in, double& E_out, double& f, int& l, int& j,
    uint64_t* seed) const;

  //! Sample an outgoing energy and evaluate the angular PDF
  //! \param[in] E_in Incoming energy in [eV]
  //! \param[in] mu Scattering cosine with respect to current direction
  //! \param[out] E_out Outgoing energy in [eV]
  //! \param[inout] seed Pseudorandom seed pointer
  //! \return Probability density for the scattering cosine
  double sample_energy_and_pdf(
    double E_in, double mu, double& E_out, uint64_t* seed) const override;

private:
  //! Secondary energy/angle distribution
  struct DistEnergySab {
    std::size_t n_e_out;              //!< Number of outgoing energies
    tensor::Tensor<double> e_out;     //!< Outgoing energies
    tensor::Tensor<double> e_out_pdf; //!< Probability density function
    tensor::Tensor<double> e_out_cdf; //!< Cumulative distribution function
    tensor::Tensor<double> mu; //!< Equiprobable angles at each outgoing energy
  };

  vector<double> energy_;              //!< Incident energies
  vector<DistEnergySab> distribution_; //!< Secondary angle-energy at
                                       //!< each incident energy
};

//==============================================================================
//! Mixed coherent/incoherent elastic angle-energy distribution
//==============================================================================

class MixedElasticAE : public AngleEnergy {
public:
  //! Construct from HDF5 file
  //
  //! \param[in] group  HDF5 group
  explicit MixedElasticAE(
    hid_t group, const CoherentElasticXS& coh_xs, const Function1D& incoh_xs);

  //! Sample distribution for an angle and energy
  //! \param[in] E_in Incoming energy in [eV]
  //! \param[out] E_out Outgoing energy in [eV]
  //! \param[out] mu Outgoing cosine with respect to current direction
  //! \param[inout] seed Pseudorandom number seed pointer
  void sample(
    double E_in, double& E_out, double& mu, uint64_t* seed) const override;

  //! Select the coherent or incoherent elastic distribution to sample
  //! \param[in] E_in Incoming energy in [eV]
  //! \param[inout] seed Pseudorandom seed pointer
  //! \return Reference to the selected angle-energy distribution
  const AngleEnergy& sample_dist(double E_in, uint64_t* seed) const;

  //! Sample an outgoing energy and evaluate the angular PDF
  //! \param[in] E_in Incoming energy in [eV]
  //! \param[in] mu Scattering cosine with respect to current direction
  //! \param[out] E_out Outgoing energy in [eV]
  //! \param[inout] seed Pseudorandom seed pointer
  //! \return Probability density for the scattering cosine
  double sample_energy_and_pdf(
    double E_in, double mu, double& E_out, uint64_t* seed) const override;

private:
  CoherentElasticAE coherent_dist_;         //!< Coherent distribution
  unique_ptr<AngleEnergy> incoherent_dist_; //!< Incoherent distribution

  const CoherentElasticXS& coherent_xs_; //!< Ref. to coherent XS
  const Function1D& incoherent_xs_;      //!< Polymorphic ref. to incoherent XS
};

//! Internal helper for evaluating a piecewise-constant PDF on discrete points.
//!
//! The underlying discrete points are represented implicitly through a
//! monotonically increasing `center(i)` function and corresponding per-point
//! `weight(i)` values. Each point contributes a rectangular bin whose
//! half-width is half the distance to its nearest neighboring center.
//!
//! \tparam CenterFn Callable returning the location of the i-th discrete value
//! \tparam WeightFn Callable returning the weight of the i-th discrete value
//! \param[in] n Number of discrete values
//! \param[in] mu_0 Point at which to evaluate the PDF
//! \param[in] a Lower bound of the domain (default: -1)
//! \param[in] b Upper bound of the domain (default: 1)
//! \return Probability density at mu_0
template<typename CenterFn, typename WeightFn>
double get_pdf_discrete_impl(std::size_t n, double mu_0, double a, double b,
  CenterFn center, WeightFn weight)
{
  if (n == 0 || mu_0 < a || mu_0 > b)
    return 0.0;

  auto evaluate_bin = [&](std::size_t i) {
    double x = center(i);
    double left_span = (i == 0) ? 2.0 * (x - a) : x - center(i - 1);
    double right_span = (i + 1 == n) ? 2.0 * (b - x) : center(i + 1) - x;
    double delta = 0.5 * std::min(left_span, right_span);
    if (delta <= 0.0)
      return 0.0;

    double left = x - delta;
    double right = x + delta;
    bool in_bin =
      (mu_0 >= left) && ((i + 1 == n) ? (mu_0 <= right) : (mu_0 < right));
    return in_bin ? weight(i) / (2.0 * delta) : 0.0;
  };

  // This is effectively a lower_bound over the sequence center(i), but the
  // sequence is implicit rather than stored in a container, so the STL
  // algorithms can not be used.
  std::size_t low = 0;
  std::size_t high = n;
  while (low < high) {
    std::size_t mid = low + (high - low) / 2;
    if (center(mid) < mu_0) {
      low = mid + 1;
    } else {
      high = mid;
    }
  }

  if (low < n) {
    double pdf = evaluate_bin(low);
    if (pdf > 0.0)
      return pdf;
  }
  if (low > 0)
    return evaluate_bin(low - 1);
  return 0.0;
}

//! Evaluate the PDF of a weighted discrete distribution at a given point.
//!
//! Given a set of discrete values mu[i] with weights w[i], this function
//! computes the probability density at mu_0 by treating each discrete value
//! as a rectangular bin. The bin half-width around each discrete value is
//! half the distance to its nearest neighbor.
//!
//! \tparam T1 Container type for discrete cosine values (must support
//!         operator[], size())
//! \tparam T2 Container type for weights (must support operator[])
//! \param[in] mu Sorted array of discrete cosine values
//! \param[in] w Weights for each discrete value (need not be normalized)
//! \param[in] mu_0 Point at which to evaluate the PDF
//! \param[in] a Lower bound of the domain (default: -1)
//! \param[in] b Upper bound of the domain (default: 1)
//! \return Probability density at mu_0
template<typename T1, typename T2>
double get_pdf_discrete(
  const T1 mu, const T2& w, double mu_0, double a = -1.0, double b = 1.0)
{
  // Returns the location of the discrete value for a given index
  auto center = [&](std::size_t i) { return mu[i]; };
  auto weight = [&](std::size_t i) { return w[i]; };
  return get_pdf_discrete_impl(mu.size(), mu_0, a, b, center, weight);
}

//! Evaluate the PDF of a discrete distribution with uniform weights
//!
//! \tparam T1 Container type for discrete cosine values
//! \param[in] mu Sorted array of discrete cosine values
//! \param[in] mu_0 Point at which to evaluate the PDF
//! \param[in] a Lower bound of the domain (default: -1)
//! \param[in] b Upper bound of the domain (default: 1)
//! \return Probability density at mu_0
template<typename T1>
double get_pdf_discrete(
  const T1 mu, double mu_0, double a = -1.0, double b = 1.0)
{
  auto center = [&](std::size_t i) { return mu[i]; };
  auto weight = [&](std::size_t i) { return 1.0 / mu.size(); };
  return get_pdf_discrete_impl(mu.size(), mu_0, a, b, center, weight);
}

//! Evaluate the PDF of a uniformly weighted distribution on interpolated points
//!
//! \tparam T1 Container type for the lower tabulated cosine values
//! \tparam T2 Container type for the upper tabulated cosine values
//! \param[in] mu0 Sorted array of discrete cosine values at the lower grid
//! \param[in] mu1 Sorted array of discrete cosine values at the upper grid
//! \param[in] f Interpolation factor between mu0 and mu1
//! \param[in] mu_0 Point at which to evaluate the PDF
//! \param[in] a Lower bound of the domain (default: -1)
//! \param[in] b Upper bound of the domain (default: 1)
//! \return Probability density at mu_0
template<typename T1, typename T2>
double get_pdf_discrete_interpolated(const T1 mu0, const T2 mu1, double f,
  double mu_0, double a = -1.0, double b = 1.0)
{
  if (mu0.size() != mu1.size())
    return 0.0;

  // Returns interpolated discrete value for a given index
  auto center = [&](std::size_t i) { return mu0[i] + f * (mu1[i] - mu0[i]); };
  auto weight = [&](std::size_t i) { return 1.0 / mu0.size(); };
  return get_pdf_discrete_impl(mu0.size(), mu_0, a, b, center, weight);
}

} // namespace openmc

#endif // OPENMC_SECONDARY_THERMAL_H
