//! \file secondary_nbody.h
//! N-body phase space distribution

#ifndef OPENMC_SECONDARY_NBODY_H
#define OPENMC_SECONDARY_NBODY_H

#include "hdf5.h"

#include "openmc/angle_energy.h"

namespace openmc {

//==============================================================================
//! Angle-energy distribution for particles emitted from neutron and
//! charged-particle reactions. This corresponds to ACE law 66 and ENDF File 6,
//! LAW=6.
//==============================================================================

class NBodyPhaseSpace : public AngleEnergy {
public:
  explicit NBodyPhaseSpace(hid_t group);

  //! Sample distribution for an angle and energy
  //! \param[in] E_in Incoming energy in [eV]
  //! \param[out] E_out Outgoing energy in [eV]
  //! \param[out] mu Outgoing cosine with respect to current direction
  //! \param[inout] seed Pseudorandom seed pointer
  void sample(double E_in, double& E_out, double& mu,
    uint64_t* seed) const override;
private:
  int n_bodies_; //!< Number of particles distributed
  double mass_ratio_; //!< Total mass of particles [neutron mass]
  double A_; //!< Atomic weight ratio
  double Q_; //!< Reaction Q-value [eV]
};

} // namespace openmc

#endif // OPENMC_SECONDARY_NBODY_H
