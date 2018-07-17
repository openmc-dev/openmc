#ifndef OPENMC_SOURCE_H
#define OPENMC_SOURCE_H

#include <memory>
#include <vector>

#include "pugixml.hpp"

#include "openmc/distribution_multi.h"
#include "openmc/distribution_spatial.h"
#include "openmc/capi.h"
#include "openmc/particle.h"

namespace openmc {

//==============================================================================
//! External source distribution
//==============================================================================

class SourceDistribution {
public:
  SourceDistribution(UPtrSpace space, UPtrAngle angle, UPtrDist energy);
  explicit SourceDistribution(pugi::xml_node node);

  Bank sample() const;
  double strength() const { return strength_; }
private:
  ParticleType particle_ {ParticleType::neutron};
  double strength_ {1.0};
  UPtrSpace space_;
  UPtrAngle angle_;
  UPtrDist energy_;
};

//==============================================================================
// Global variables
//==============================================================================

extern std::vector<SourceDistribution> external_sources;

//==============================================================================
// Functions
//==============================================================================

//! Initialize source bank from file/distribution
extern "C" void initialize_source();

//! Sample a site from all external source distributions in proportion to their
//! source strength
extern "C" Bank sample_external_source();

} // namespace openmc

#endif // OPENMC_SOURCE_H
