#ifndef OPENMC_ANGLE_ENERGY_H
#define OPENMC_ANGLE_ENERGY_H

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
  virtual void sample(
    double E_in, double& E_out, double& mu, uint64_t* seed) const = 0;
  virtual double sample_energy_and_pdf(
    double E_in, double mu, double& E_out, uint64_t* seed) const = 0;
  virtual ~AngleEnergy() = default;
};

} // namespace openmc

#endif // OPENMC_ANGLE_ENERGY_H
