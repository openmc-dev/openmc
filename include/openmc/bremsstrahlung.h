#ifndef OPENMC_BREMSSTRAHLUNG_H
#define OPENMC_BREMSSTRAHLUNG_H

#include "openmc/particle.h"

#include <gsl/gsl>
#include "xtensor/xtensor.hpp"

namespace openmc {

//==============================================================================
// Bremsstrahlung classes
//==============================================================================

class BremsstrahlungData {
public:
  // Methods
  void copy_to_device();
  void release_from_device();

  // Data
  xt::xtensor<double, 2> pdf_; //!< Bremsstrahlung energy PDF
  xt::xtensor<double, 2> cdf_; //!< Bremsstrahlung energy CDF
  xt::xtensor<double, 1> yield_; //!< Photon yield

  double* device_pdf_ {nullptr};
  double* device_cdf_ {nullptr};
  double* device_yield_ {nullptr};
  double pdf(gsl::index i, gsl::index j) const;
  double cdf(gsl::index i, gsl::index j) const;
};

class Bremsstrahlung {
public:
  // Methods
  void copy_to_device();
  void release_from_device();

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

#pragma omp declare target
extern double* device_ttb_e_grid;
extern size_t ttb_e_grid_size;
#pragma omp end declare target

} // namespace data

//==============================================================================
// Non-member functions
//==============================================================================

#pragma omp declare target
void thick_target_bremsstrahlung(Particle& p, double* E_lost);
#pragma omp end declare target

} // namespace openmc

#endif // OPENMC_BREMSSTRAHLUNG_H
