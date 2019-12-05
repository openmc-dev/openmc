//! \file reaction_product.h
//! Data for a reaction product

#ifndef OPENMC_REACTION_PRODUCT_H
#define OPENMC_REACTION_PRODUCT_H

#include <memory> // for unique_ptr
#include <vector> // for vector

#include "hdf5.h"

#include "openmc/angle_energy.h"
#include "openmc/endf.h"
#include "openmc/particle.h"

namespace openmc {

//==============================================================================
//! Data for a reaction product including its yield and angle-energy
//! distributions, each of which has a given probability of occurring for a
//! given incoming energy. In general, most products only have one angle-energy
//! distribution, but for some cases (e.g., (n,2n) in certain nuclides) multiple
//! distinct distributions exist.
//==============================================================================

class ReactionProduct {
public:
  //! Emission mode for product
  enum class EmissionMode {
    prompt,  // Prompt emission of secondary particle
    delayed, // Yield represents total emission (prompt + delayed)
    total    // Delayed emission of secondary particle
  };

  using Secondary = std::unique_ptr<AngleEnergy>;

  //! Construct reaction product from HDF5 data
  //! \param[in] group HDF5 group containing data
  explicit ReactionProduct(hid_t group);

  //! Sample an outgoing angle and energy
  //! \param[in] E_in Incoming energy in [eV]
  //! \param[out] E_out Outgoing energy in [eV]
  //! \param[out] mu Outgoing cosine with respect to current direction
  //! \param[inout] seed Pseudorandom seed pointer
  void sample(double E_in, double& E_out, double& mu, uint64_t* seed) const;

  Particle::Type particle_; //!< Particle type
  EmissionMode emission_mode_; //!< Emission mode
  double decay_rate_; //!< Decay rate (for delayed neutron precursors) in [1/s]
  std::unique_ptr<Function1D> yield_; //!< Yield as a function of energy
  std::vector<Tabulated1D> applicability_; //!< Applicability of distribution
  std::vector<Secondary> distribution_; //!< Secondary angle-energy distribution
};

} // namespace opemc

#endif // OPENMC_REACTION_PRODUCT_H
