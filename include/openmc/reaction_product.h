//! \file reaction_product.h
//! Data for a reaction product

#ifndef OPENMC_REACTION_PRODUCT_H
#define OPENMC_REACTION_PRODUCT_H

#include <memory> // for unique_ptr
#include <vector> // for vector

#include "hdf5.h"
#include <gsl/gsl>

#include "openmc/angle_energy.h"
#include "openmc/endf.h"
#include "openmc/endf_flat.h"
#include "openmc/particle.h"
#include "openmc/secondary_flat.h"
#include "openmc/serialize.h"

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

  void serialize(DataBuffer& buffer) const;

  Particle::Type particle_; //!< Particle type
  EmissionMode emission_mode_; //!< Emission mode
  double decay_rate_; //!< Decay rate (for delayed neutron precursors) in [1/s]
  std::unique_ptr<Function1D> yield_; //!< Yield as a function of energy
  std::vector<Tabulated1D> applicability_; //!< Applicability of distribution
  std::vector<Secondary> distribution_; //!< Secondary angle-energy distribution
};

class ReactionProductFlat {
public:
  // Constructors
  #pragma omp declare target
  explicit ReactionProductFlat(const uint8_t* data);
  #pragma omp end declare target

  void sample(double E_in, double& E_out, double& mu, uint64_t* seed) const;

  Particle::Type particle() const;
  ReactionProduct::EmissionMode emission_mode() const;
  double decay_rate() const;
  Function1DFlat yield() const;
  Tabulated1DFlat applicability(gsl::index i) const;
  AngleEnergyFlat distribution(gsl::index i) const;

  size_t n_distribution() const { return n_distribution_; }

private:
  // Data members
  const uint8_t* data_;
  size_t yield_size_;
  size_t n_applicability_;
  size_t n_distribution_;
};

class ReactionProductFlatContainer {
public:
  // Constructors
  explicit ReactionProductFlatContainer(const ReactionProduct& product);

  #pragma omp declare target
  void sample(double E_in, double& E_out, double& mu, uint64_t* seed) const;
  #pragma omp end declare target

  void copy_to_device();
  void release_from_device();

  ReactionProductFlat obj() const;

  Particle::Type particle() const { return this->obj().particle(); }
  ReactionProduct::EmissionMode emission_mode() const { return this->obj().emission_mode(); }
  double decay_rate() const { return this->obj().decay_rate(); }
  Function1DFlat yield() const { return this->obj().yield(); }
  AngleEnergyFlat distribution(gsl::index i) const { return this->obj().distribution(i); }

  const uint8_t* data() const { return buffer_.data_; }
private:
  // Data members
  DataBuffer buffer_;
};

} // namespace opemc

#endif // OPENMC_REACTION_PRODUCT_H
