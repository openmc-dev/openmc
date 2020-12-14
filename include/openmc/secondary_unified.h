//! \file secondary_unified.h
//! Unified angle-energy distribution

#ifndef OPENMC_SECONDARY_UNIFIED_H
#define OPENMC_SECONDARY_UNIFIED_H

#include <cstdint> // for uint8_t
#include <memory> // for unique_ptr

#include "openmc/angle_energy.h"

namespace openmc {

enum AngleEnergyType {
  UNCORRELATED,
  KALBACH_MANN,
  CORRELATED,
  NBODY
};

class DataBuffer {
public:
  explicit DataBuffer(size_t n);

  void add(int value);
  void add(double value);

  std::unique_ptr<uint8_t[]> data_;
  size_t offset_{0};
};

class UnifiedAngleEnergy : public AngleEnergy {
public:
  // Constructors
  UnifiedAngleEnergy(AngleEnergyType type, DataBuffer buffer);

  //! Sample distribution for an angle and energy
  //! \param[in] E_in Incoming energy in [eV]
  //! \param[out] E_out Outgoing energy in [eV]
  //! \param[out] mu Outgoing cosine with respect to current direction
  //! \param[inout] seed Pseudorandom seed pointer
  void sample(double E_in, double& E_out, double& mu, uint64_t* seed) const override;

  uint8_t* data() { return buffer_.data_.get(); }
private:
  // Data members
  AngleEnergyType type_;
  DataBuffer buffer_;
};


} // namespace openmc

#endif // OPENMC_SECONDARY_UNIFIED_H
