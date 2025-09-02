#ifndef OPENMC_ANGLE_ENERGY_H
#define OPENMC_ANGLE_ENERGY_H

#include <cstdint>
#include "openmc/error.h"

namespace openmc {

//==============================================================================
//! Abstract type that defines a correlated or uncorrelated angle-energy
//! distribution that is a function of incoming energy. Each derived type must
//! implement a sample() method that returns an outgoing energy and
//! scattering cosine given an incoming energy.
//==============================================================================

class AngleEnergy {
public:
  virtual void sample(
    double E_in, double& E_out, double& mu, uint64_t* seed) const = 0;
  virtual double get_pdf(double, double&, double, uint64_t*) const
  {
  fatal_error("get_pdf not available for this AngleEnergy type");
  }
  virtual ~AngleEnergy() = default;
};

} // namespace openmc

#endif // OPENMC_ANGLE_ENERGY_H
