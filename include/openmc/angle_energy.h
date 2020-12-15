#ifndef OPENMC_ANGLE_ENERGY_H
#define OPENMC_ANGLE_ENERGY_H

#include "openmc/random_lcg.h"
#include <cstdint>

namespace openmc {

//==============================================================================
//! Abstract type that defines a correlated or uncorrelated angle-energy
//! distribution that is a function of incoming energy. Each derived type must
//! implement a sample() method that returns an outgoing energy and
//! scattering cosine given an incoming energy.
//==============================================================================

class AngleEnergy {
public:
  HD virtual void sample(xsfloat E_in, xsfloat& E_out, xsfloat& mu,
    uint64_t* seed) const = 0;
  AngleEnergy(AngleEnergy&&) = default;
  AngleEnergy() = default;
  virtual ~AngleEnergy() = default;

  // Center of mass mu sampling
  HD virtual double sample_mu(double const& E, uint64_t* seed)
  {
    return 2.0 * prn(seed) - 1.0;
  }
};

}

#endif // OPENMC_ANGLE_ENERGY_H
