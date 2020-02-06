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
private:
  const CoherentElasticXS& xs_; //!< Coherent elastic scattering cross section
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
  explicit IncoherentElasticAEDiscrete(hid_t group, const std::vector<double>& energy);

  //! Sample distribution for an angle and energy
  //! \param[in] E_in Incoming energy in [eV]
  //! \param[out] E_out Outgoing energy in [eV]
  //! \param[out] mu Outgoing cosine with respect to current direction
  //! \param[inout] seed Pseudorandom number seed pointer
  void sample(double E_in, double& E_out, double& mu,
    uint64_t* seed) const override;
private:
  const std::vector<double>& energy_; //!< Energies at which cosines are tabulated
  xt::xtensor<double, 2> mu_out_; //!< Cosines for each incident energy
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
private:
  const std::vector<double>& energy_; //!< Incident energies
  xt::xtensor<double, 2> energy_out_; //!< Outgoing energies for each incident energy
  xt::xtensor<double, 3> mu_out_; //!< Outgoing cosines for each incident/outgoing energy
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
  void sample(double E_in, double& E_out, double& mu,
    uint64_t* seed) const override;
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


} // namespace openmc

#endif // OPENMC_SECONDARY_THERMAL_H
