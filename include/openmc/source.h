//! \file source.h
//! \brief External source distributions

#ifndef OPENMC_SOURCE_H
#define OPENMC_SOURCE_H

#include <memory>
#include <vector>

#include "pugixml.hpp"

#include "openmc/distribution_multi.h"
#include "openmc/distribution_spatial.h"
#include "openmc/particle.h"

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

class SourceDistribution;

namespace model {

extern std::vector<SourceDistribution> external_sources;

} // namespace model

//==============================================================================
//! External source distribution
//==============================================================================

class SourceDistribution {
public:
  // Constructors
  SourceDistribution(UPtrSpace space, UPtrAngle angle, UPtrDist energy);
  explicit SourceDistribution(pugi::xml_node node);

  //! Sample from the external source distribution
  //! \param[inout] seed Pseudorandom seed pointer
  //! \return Sampled site
  Particle::Bank sample(uint64_t* seed) const;

  // Properties
  Particle::Type particle_type() const { return particle_; }
  double strength() const { return strength_; }

  // Make observing pointers available
  SpatialDistribution* space() const { return space_.get(); }
  UnitSphereDistribution* angle() const { return angle_.get(); }
  Distribution* energy() const { return energy_.get(); }

private:
  Particle::Type particle_ {Particle::Type::neutron}; //!< Type of particle emitted
  double strength_ {1.0}; //!< Source strength
  UPtrSpace space_; //!< Spatial distribution
  UPtrAngle angle_; //!< Angular distribution
  UPtrDist energy_; //!< Energy distribution
};

//==============================================================================
// Functions
//==============================================================================

//! Initialize source bank from file/distribution
extern "C" void initialize_source();

//! Sample a site from all external source distributions in proportion to their
//! source strength
//! \param[inout] seed Pseudorandom seed pointer
//! \return Sampled source site
Particle::Bank sample_external_source(uint64_t* seed);

//! Sample a site from custom source library
Particle::Bank sample_custom_source_library(uint64_t* seed);

//! Load custom source library
void load_custom_source_library();

//! Release custom source library
void close_custom_source_library();

//! Fill source bank at the end of a generation for dlopen based source simulation
void fill_source_bank_custom_source();

void free_memory_source();

} // namespace openmc

#endif // OPENMC_SOURCE_H
