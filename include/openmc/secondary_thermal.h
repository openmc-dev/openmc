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

//! Helper for returning a uniform weight via operator[] regardless of index
struct DoubleVector {
  double data;
  const double& operator[](size_t index) const { return data; }
};

//! Evaluate the PDF of a discrete distribution at a given point.
//!
//! Given a set of discrete values mu[i] with weights w[i], this function
//! computes the probability density at mu_0 by treating each discrete value
//! as a rectangular bin. The bin half-width around each discrete value is
//! half the distance to its nearest neighbor.
//!
//! \tparam T1 Container type for discrete cosine values (must support
//!         operator[], begin(), end(), size())
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
  // Clamp mu_0 to the domain [a, b]
  if (mu_0 < a)
    mu_0 = a;
  if (mu_0 > b)
    mu_0 = b;

  // Find the nearest discrete values bracketing mu_0:
  //   a0 = closest discrete value <= mu_0 (or domain bound a)
  //   a1 = the discrete value before a0 (or domain bound a)
  //   b0 = closest discrete value >= mu_0 (or domain bound b)
  //   b1 = the discrete value after b0 (or domain bound b)
  double a0, a1, b0, b1;
  int32_t ai = -1;
  int32_t bi = -1;
  if (mu_0 > mu[0]) {
    ai = lower_bound_index(mu.begin(), mu.end(), mu_0);
    a0 = mu[ai];
    a1 = (ai > 1) ? mu[ai - 1] : a;
  } else {
    a0 = a;
    a1 = a;
  }
  if (mu_0 < mu[mu.size() - 1]) {
    bi = upper_bound_index(mu.begin(), mu.end(), mu_0);
    b0 = mu[bi];
    b1 = (bi < mu.size() - 1) ? mu[bi + 1] : b;
  } else {
    b0 = b;
    b1 = b;
  }

  // Compute the half-width of the rectangular bin around each bracketing
  // discrete value. The half-width is half the distance to its nearest
  // neighbor (or domain boundary).
  double delta_a = 0.5 * std::min(b0 - a0, a0 - a1);
  double delta_b = 0.5 * std::min(b1 - b0, b0 - a0);

  // Determine which bin mu_0 falls in and return the corresponding density
  if (mu_0 < a0 + delta_a)
    return w[ai] / (2.0 * delta_a);
  else if (mu_0 + delta_b < b0)
    return w[bi] / (2.0 * delta_b);
  else
    return 0.0;
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
  DoubleVector w {1.0 / mu.size()};
  return get_pdf_discrete(mu, w, mu_0, a, b);
}

} // namespace openmc

#endif // OPENMC_SECONDARY_THERMAL_H
