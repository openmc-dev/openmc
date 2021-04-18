#ifndef OPENMC_ANGLE_ENERGY_H
#define OPENMC_ANGLE_ENERGY_H

#include <cstdint>

#include "openmc/serialize.h"

namespace openmc {

//==============================================================================
//! Abstract type that defines a correlated or uncorrelated angle-energy
//! distribution that is a function of incoming energy. Each derived type must
//! implement a sample() method that returns an outgoing energy and
//! scattering cosine given an incoming energy.
//==============================================================================

class AngleEnergy {
public:
  virtual void sample(double E_in, double& E_out, double& mu,
    uint64_t* seed) const = 0;
  virtual ~AngleEnergy() = default;

  virtual void serialize(DataBuffer& buffer) const = 0;
};

}

#endif // OPENMC_ANGLE_ENERGY_H
