//! \file secondary_thermal.h
//! Angle-energy distributions for thermal scattering

#ifndef OPENMC_SECONDARY_THERMAL_H
#define OPENMC_SECONDARY_THERMAL_H

#include "openmc/angle_energy.h"
#include "openmc/endf.h"
#include "openmc/secondary_correlated.h"
#include "openmc/tensor.h"
#include "openmc/vector.h"

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
  HD void sample(xsfloat E_in, xsfloat& E_out, xsfloat& mu,
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
  HD void sample(xsfloat E_in, xsfloat& E_out, xsfloat& mu,
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
  explicit IncoherentElasticAEDiscrete(
    hid_t group, const vector<xsfloat>& energy);

  //! Sample distribution for an angle and energy
  //! \param[in] E_in Incoming energy in [eV]
  //! \param[out] E_out Outgoing energy in [eV]
  //! \param[out] mu Outgoing cosine with respect to current direction
  //! \param[inout] seed Pseudorandom number seed pointer
  HD void sample(xsfloat E_in, xsfloat& E_out, xsfloat& mu,
    uint64_t* seed) const override;
private:
  const vector<xsfloat>& energy_;  //!< Energies at which cosines are tabulated
  tensor<xsfloat, 2> mu_out_;      //!< Cosines for each incident energy
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
    hid_t group, const vector<xsfloat>& energy);

  //! Sample distribution for an angle and energy
  //! \param[in] E_in Incoming energy in [eV]
  //! \param[out] E_out Outgoing energy in [eV]
  //! \param[out] mu Outgoing cosine with respect to current direction
  //! \param[inout] seed Pseudorandom number seed pointer
  HD void sample(xsfloat E_in, xsfloat& E_out, xsfloat& mu,
    uint64_t* seed) const override;
private:
  const vector<xsfloat>& energy_;      //!< Incident energies
  tensor<xsfloat, 2>
    energy_out_; //!< Outgoing energies for each incident energy
  tensor<xsfloat, 3>
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
  HD void sample(xsfloat E_in, xsfloat& E_out, xsfloat& mu,
    uint64_t* seed) const override;
private:
  //! Secondary energy/angle distribution
  struct DistEnergySab {
    std::size_t n_e_out; //!< Number of outgoing energies
    vector<xsfloat> e_out;     //!< Outgoing energies
    vector<xsfloat> e_out_pdf; //!< Probability density function
    vector<xsfloat> e_out_cdf; //!< Cumulative distribution function
    tensor<xsfloat, 2> mu;     //!< Equiprobable angles at each outgoing energy
  };

  vector<xsfloat> energy_;              //!< Incident energies
  vector<DistEnergySab> distribution_; //!< Secondary angle-energy at
                                       //!< each incident energy
};


} // namespace openmc

#endif // OPENMC_SECONDARY_THERMAL_H
