//! \file secondary_flat.h
//! Unified angle-energy distribution

#ifndef OPENMC_SECONDARY_FLAT_H
#define OPENMC_SECONDARY_FLAT_H

#include "openmc/angle_energy.h"
#include "openmc/serialize.h"

namespace openmc {

enum class AngleEnergyType {
  UNCORRELATED,
  KALBACH_MANN,
  CORRELATED,
  NBODY
};

class AngleEnergyFlat {
public:
  // Constructors
  explicit AngleEnergyFlat(const AngleEnergy& dist);
  explicit AngleEnergyFlat(DataBuffer buffer);

  //! Sample distribution for an angle and energy
  //! \param[in] E_in Incoming energy in [eV]
  //! \param[out] E_out Outgoing energy in [eV]
  //! \param[out] mu Outgoing cosine with respect to current direction
  //! \param[inout] seed Pseudorandom seed pointer
  void sample(double E_in, double& E_out, double& mu, uint64_t* seed) const;

  const uint8_t* data() const { return buffer_.data_.get(); }
  AngleEnergyType type() const { return type_; }
private:
  // Data members
  AngleEnergyType type_;
  DataBuffer buffer_;
};

} // namespace openmc

#endif // OPENMC_SECONDARY_FLAT_H
