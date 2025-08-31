//! \file reaction_product.h
//! Data for a reaction product

#ifndef OPENMC_REACTION_PRODUCT_H
#define OPENMC_REACTION_PRODUCT_H

#include "hdf5.h"

#include "openmc/angle_energy.h"
#include "openmc/chain.h"
#include "openmc/endf.h"
#include "openmc/memory.h" // for unique_ptr
#include "openmc/particle.h"
#include "openmc/tallies/filter.h"
#include "openmc/tallies/tally.h"
#include "openmc/vector.h" // for vector

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

  using Secondary = unique_ptr<AngleEnergy>;

  //! Construct reaction product from HDF5 data
  //! \param[in] group HDF5 group containing data
  explicit ReactionProduct(hid_t group);

  //! Construct reaction product for decay photon from chain nuclide product
  //! \param[in] product Chain nuclide product
  explicit ReactionProduct(const ChainNuclide::Product& product);

  //! Sample an outgoing angle and energy
  //! \param[in] E_in Incoming energy in [eV]
  //! \param[out] E_out Outgoing energy in [eV]
  //! \param[out] mu Outgoing cosine with respect to current direction
  //! \param[inout] seed Pseudorandom seed pointer
  void sample(double E_in, double& E_out, double& mu, uint64_t* seed) const;

  void get_pdf(int i_tally, double E_in, double& E_out, uint64_t* seed,
    Particle& p, std::vector<double>& mu_cm, std::vector<double>& Js,
    std::vector<Particle>& ghost_particles,
    std::vector<double>& pdfs_lab) const;

  ParticleType particle_;      //!< Particle type
  EmissionMode emission_mode_; //!< Emission mode
  double decay_rate_; //!< Decay rate (for delayed neutron precursors) in [1/s]
  unique_ptr<Function1D> yield_;      //!< Yield as a function of energy
  vector<Tabulated1D> applicability_; //!< Applicability of distribution
  vector<Secondary> distribution_;    //!< Secondary angle-energy distribution
  int parent_nuclide_ = -1;           //!< Index of chain nuclide that is parent
};

} // namespace openmc

#endif // OPENMC_REACTION_PRODUCT_H
