#ifndef OPENMC_BREMSSTRAHLUNG_H
#define OPENMC_BREMSSTRAHLUNG_H

#include "openmc/particle.h"

#include "openmc/tensor.h"

namespace openmc {

//==============================================================================
// Bremsstrahlung classes
//==============================================================================

class BremsstrahlungData {
public:
  // Data
  tensor::Tensor<double> pdf;   //!< Bremsstrahlung energy PDF
  tensor::Tensor<double> cdf;   //!< Bremsstrahlung energy CDF
  tensor::Tensor<double> yield; //!< Photon yield
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

extern tensor::Tensor<double>
  ttb_e_grid; //! energy T of incident electron in [eV]
extern tensor::Tensor<double>
  ttb_k_grid; //! reduced energy W/T of emitted photon

} // namespace data

//==============================================================================
// Global variables
//==============================================================================

void thick_target_bremsstrahlung(Particle& p, double* E_lost);

} // namespace openmc

#endif // OPENMC_BREMSSTRAHLUNG_H
