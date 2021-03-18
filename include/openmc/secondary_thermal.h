//! \file secondary_thermal.h
//! Angle-energy distributions for thermal scattering

#ifndef OPENMC_SECONDARY_THERMAL_H
#define OPENMC_SECONDARY_THERMAL_H

#include "openmc/angle_energy.h"
#include "openmc/endf.h"
#include "openmc/secondary_correlated.h"

#include <hdf5.h>
#include "xtensor/xtensor.hpp"

#include <vector>

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
  void sample(double E_in, double& E_out, double& mu,
    uint64_t* seed) const override;

  void serialize(DataBuffer& buffer) const override;
private:
  const CoherentElasticXS& xs_; //!< Coherent elastic scattering cross section
};

class CoherentElasticAEFlat {
public:
  #pragma omp declare target
  explicit CoherentElasticAEFlat(const uint8_t* data) : data_(data) { }

  //! Sample distribution for an angle and energy
  //! \param[in] E_in Incoming energy in [eV]
  //! \param[out] E_out Outgoing energy in [eV]
  //! \param[out] mu Outgoing cosine with respect to current direction
  //! \param[inout] seed Pseudorandom seed pointer
  void sample(double E_in, double& E_out, double& mu, uint64_t* seed) const;
  #pragma omp end declare target

private:
  CoherentElasticXSFlat xs() const;

  const uint8_t* data_;
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
  void sample(double E_in, double& E_out, double& mu,
    uint64_t* seed) const override;

  void serialize(DataBuffer& buffer) const override;
private:
  double debye_waller_;
};

class IncoherentElasticAEFlat {
public:
  #pragma omp declare target
  explicit IncoherentElasticAEFlat(const uint8_t* data) : data_(data) { }

  //! Sample distribution for an angle and energy
  //! \param[in] E_in Incoming energy in [eV]
  //! \param[out] E_out Outgoing energy in [eV]
  //! \param[out] mu Outgoing cosine with respect to current direction
  //! \param[inout] seed Pseudorandom number seed pointer
  void sample(double E_in, double& E_out, double& mu, uint64_t* seed) const;
  #pragma omp end declare target

private:
  double debye_waller() const;

  const uint8_t* data_;
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
  explicit IncoherentElasticAEDiscrete(hid_t group, const std::vector<double>& energy);

  //! Sample distribution for an angle and energy
  //! \param[in] E_in Incoming energy in [eV]
  //! \param[out] E_out Outgoing energy in [eV]
  //! \param[out] mu Outgoing cosine with respect to current direction
  //! \param[inout] seed Pseudorandom number seed pointer
  void sample(double E_in, double& E_out, double& mu,
    uint64_t* seed) const override;

  void serialize(DataBuffer& buffer) const override;
private:
  const std::vector<double>& energy_; //!< Energies at which cosines are tabulated
  xt::xtensor<double, 2> mu_out_; //!< Cosines for each incident energy
};

class IncoherentElasticAEDiscreteFlat {
public:
  #pragma omp declare target
  explicit IncoherentElasticAEDiscreteFlat(const uint8_t* data);

  //! Sample distribution for an angle and energy
  //! \param[in] E_in Incoming energy in [eV]
  //! \param[out] E_out Outgoing energy in [eV]
  //! \param[out] mu Outgoing cosine with respect to current direction
  //! \param[inout] seed Pseudorandom number seed pointer
  void sample(double E_in, double& E_out, double& mu, uint64_t* seed) const;
  #pragma omp end declare target

private:
  gsl::span<const double> energy() const;
  double mu_out(gsl::index i, gsl::index j) const;

  const uint8_t* data_;
  size_t n_e_out_;
  size_t n_mu_;
  const double* mu_out_;
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
  explicit IncoherentInelasticAEDiscrete(hid_t group, const std::vector<double>& energy);

  //! Sample distribution for an angle and energy
  //! \param[in] E_in Incoming energy in [eV]
  //! \param[out] E_out Outgoing energy in [eV]
  //! \param[out] mu Outgoing cosine with respect to current direction
  //! \param[inout] seed Pseudorandom number seed pointer
  void sample(double E_in, double& E_out, double& mu,
    uint64_t* seed) const override;

  void serialize(DataBuffer& buffer) const override;
private:
  const std::vector<double>& energy_; //!< Incident energies
  xt::xtensor<double, 2> energy_out_; //!< Outgoing energies for each incident energy
  xt::xtensor<double, 3> mu_out_; //!< Outgoing cosines for each incident/outgoing energy
  bool skewed_; //!< Whether outgoing energy distribution is skewed
};

class IncoherentInelasticAEDiscreteFlat {
public:
  #pragma omp declare target
  explicit IncoherentInelasticAEDiscreteFlat(const uint8_t* data);

  //! Sample distribution for an angle and energy
  //! \param[in] E_in Incoming energy in [eV]
  //! \param[out] E_out Outgoing energy in [eV]
  //! \param[out] mu Outgoing cosine with respect to current direction
  //! \param[inout] seed Pseudorandom number seed pointer
  void sample(double E_in, double& E_out, double& mu, uint64_t* seed) const;
  #pragma omp end declare target

private:
  gsl::span<const double> energy() const;
  double energy_out(gsl::index i, gsl::index j) const;
  double mu_out(gsl::index i, gsl::index j, gsl::index k) const;
  bool skewed() const;

  const uint8_t* data_;
  size_t n_energy_;
  size_t n_e_out_;
  size_t n_mu_;
  const double* energy_out_;
  const double* mu_out_;
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
  void sample(double E_in, double& E_out, double& mu,
    uint64_t* seed) const override;

  void serialize(DataBuffer& buffer) const override;
private:
  //! Secondary energy/angle distribution
  struct DistEnergySab {
    std::size_t n_e_out; //!< Number of outgoing energies
    xt::xtensor<double, 1> e_out;     //!< Outgoing energies
    xt::xtensor<double, 1> e_out_pdf; //!< Probability density function
    xt::xtensor<double, 1> e_out_cdf; //!< Cumulative distribution function
    xt::xtensor<double, 2> mu; //!< Equiprobable angles at each outgoing energy
  };

  std::vector<double> energy_; //!< Incident energies
  std::vector<DistEnergySab> distribution_; //!< Secondary angle-energy at
                                            //!< each incident energy
};

class DistEnergySabFlat {
public:
  #pragma omp declare target
  explicit DistEnergySabFlat(const uint8_t* data);
  #pragma omp end declare target

  size_t n_e_out() const { return n_e_out_; }
  size_t n_mu() const { return n_mu_; }
  gsl::span<const double> e_out() const;
  gsl::span<const double> e_out_pdf() const;
  gsl::span<const double> e_out_cdf() const;
  double mu(gsl::index i, gsl::index j) const;
private:
  const uint8_t* data_;
  size_t n_e_out_;
  size_t n_mu_;
};

class IncoherentInelasticAEFlat{
public:
  #pragma omp declare target
  explicit IncoherentInelasticAEFlat(const uint8_t* data);

  //! Sample distribution for an angle and energy
  //! \param[in] E_in Incoming energy in [eV]
  //! \param[out] E_out Outgoing energy in [eV]
  //! \param[out] mu Outgoing cosine with respect to current direction
  //! \param[inout] seed Pseudorandom number seed pointer
  void sample(double E_in, double& E_out, double& mu, uint64_t* seed) const;
  #pragma omp end declare target
private:
  gsl::span<const double> energy() const;
  DistEnergySabFlat distribution(gsl::index i) const;

  const uint8_t* data_;
  size_t n_energy_;
};


} // namespace openmc

#endif // OPENMC_SECONDARY_THERMAL_H
