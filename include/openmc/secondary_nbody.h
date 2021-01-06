//! \file secondary_nbody.h
//! N-body phase space distribution

#ifndef OPENMC_SECONDARY_NBODY_H
#define OPENMC_SECONDARY_NBODY_H

#include <cstdint>

#include <hdf5.h>

#include "openmc/angle_energy.h"
#include "openmc/secondary_unified.h"

namespace openmc {

//==============================================================================
//! Angle-energy distribution for particles emitted from neutron and
//! charged-particle reactions. This corresponds to ACE law 66 and ENDF File 6,
//! LAW=6.
//==============================================================================

class NBodyPhaseSpace : public AngleEnergy {
public:
  explicit NBodyPhaseSpace(hid_t group);

  NBodyPhaseSpace(int n_bodies, double mass_ratio, double A, double Q)
    : n_bodies_(n_bodies), mass_ratio_(mass_ratio), A_(A), Q_(Q) { }

  //! Sample distribution for an angle and energy
  //! \param[in] E_in Incoming energy in [eV]
  //! \param[out] E_out Outgoing energy in [eV]
  //! \param[out] mu Outgoing cosine with respect to current direction
  //! \param[inout] seed Pseudorandom seed pointer
  void sample(double E_in, double& E_out, double& mu, uint64_t* seed) const override;

  void serialize(DataBuffer& buffer) const override;
private:
  int n_bodies_; //!< Number of particles distributed
  double mass_ratio_; //!< Total mass of particles [neutron mass]
  double A_; //!< Atomic weight ratio
  double Q_; //!< Reaction Q-value [eV]
};

class NBodyPhaseSpaceFlat {
public:
  explicit NBodyPhaseSpaceFlat(const uint8_t* data) : data_(data) { }

  void sample(double E_in, double& E_out, double& mu, uint64_t* seed) const;

  int n_bodies() const { return *reinterpret_cast<const int*>(data_ + 4); };
  double mass_ratio() const { return *reinterpret_cast<const double*>(data_ + 8); }
  double A() const { return *reinterpret_cast<const double*>(data_ + 16); }
  double Q() const { return *reinterpret_cast<const double*>(data_ + 24); }
private:
  const uint8_t* data_;
};

} // namespace openmc

#endif // OPENMC_SECONDARY_NBODY_H
