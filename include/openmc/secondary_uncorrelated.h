//! \file secondary_uncorrelated.h
//! Uncorrelated angle-energy distribution

#ifndef OPENMC_SECONDARY_UNCORRELATED_H
#define OPENMC_SECONDARY_UNCORRELATED_H

#include "hdf5.h"

#include "openmc/angle_energy.h"
#include "openmc/distribution_angle.h"
#include "openmc/distribution_energy.h"
#include "openmc/memory.h"
#include "openmc/vector.h"

namespace openmc {

//==============================================================================
//! Uncorrelated angle-energy distribution. This corresponds to when an energy
//! distribution is given in ENDF File 5/6 and an angular distribution is given
//! in ENDF File 4.
//==============================================================================

class UncorrelatedAngleEnergy : public AngleEnergy {
public:
  explicit UncorrelatedAngleEnergy(hid_t group);

  //! Sample distribution for an angle and energy
  //! \param[in] E_in Incoming energy in [eV]
  //! \param[out] E_out Outgoing energy in [eV]
  //! \param[out] mu Outgoing cosine with respect to current direction
  //! \param[inout] seed Pseudorandom seed pointer
  HD void sample(xsfloat E_in, xsfloat& E_out, xsfloat& mu,
    uint64_t* seed) const override;

  HD virtual double sample_mu(double const& E, uint64_t* seed) override
  {
    return angle_.sample(E, seed);
  }

  // Accessors
  AngleDistribution& angle() { return angle_; }
private:
  AngleDistribution angle_; //!< Angle distribution
  unique_ptr<EnergyDistribution> energy_; //!< Energy distribution
};

} // namespace openmc

#endif // OPENMC_SECONDARY_UNCORRELATED_H
