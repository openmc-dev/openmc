//! \file source.h
//! \brief External source distributions

#ifndef OPENMC_SOURCE_H
#define OPENMC_SOURCE_H

#include <memory>
#include <vector>
#ifdef DAGMC
#include "pyne_source_sampling/source_sampling.h"
#endif

#include "pugixml.hpp"

#include "openmc/distribution_multi.h"
#include "openmc/distribution_spatial.h"
#include "openmc/particle.h"
#include "pyne/source_sampling.h"

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
  //! \return Sampled site
  Particle::Bank sample() const;

  // Properties
  double strength() const { return strength_; }
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
//! \return Sampled source site
Particle::Bank sample_external_source();

#ifdef DAGMC
//! Initialize pyne sampler instance
//! \return Sampler*
pyne::Sampler* initialize_pyne_sampler();

//! Sample a site from pyne source
//! \return Sampled source site
Particle::Bank sample_pyne_source(pyne::Sampler*);

//! Convert a pyne source particle to an openmc source site
//! \return Sampled source site
Particle::Bank convert_pyne_source_particle(pyne::SourceParticle);

//! Check the source particle converted from pyne source
//! \return Bool
bool check_pyne_source_particle(Particle::Bank);
#endif

//! Fill source bank at end of generation for fixed source simulations
void fill_source_bank_fixedsource();

void free_memory_source();

} // namespace openmc

#endif // OPENMC_SOURCE_H
