#ifndef OPENMC_BREMSSTRAHLUNG_H
#define OPENMC_BREMSSTRAHLUNG_H

#include "openmc/particle.h"

#include "xtensor/xtensor.hpp"

namespace openmc {

//==============================================================================
// Bremsstrahlung classes
//==============================================================================

class BremsstrahlungData {
public:
  // Data
  xt::xtensor<double, 2> pdf; //!< Bremsstrahlung energy PDF
  xt::xtensor<double, 2> cdf; //!< Bremsstrahlung energy CDF
  xt::xtensor<double, 1> yield; //!< Photon yield
};

class Bremsstrahlung {
public:
  // Data
  BremsstrahlungData electron;
  BremsstrahlungData positron;
};

//==============================================================================
// Global variables
//==============================================================================

namespace data {

extern xt::xtensor<double, 1> ttb_e_grid; //! energy T of incident electron in [eV]
extern xt::xtensor<double, 1> ttb_k_grid; //! reduced energy W/T of emitted photon

} // namespace data

//==============================================================================
// Global variables
//==============================================================================

void thick_target_bremsstrahlung(Particle& p, double* E_lost);

} // namespace openmc

#endif // OPENMC_BREMSSTRAHLUNG_H
