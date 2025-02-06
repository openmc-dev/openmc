#ifndef OPENMC_IFP_H
#define OPENMC_IFP_H

#include "openmc/particle.h"
#include "openmc/particle_data.h"

namespace openmc {

//! Iterated Fission Probability (IFP) method
void ifp(Particle& p, SourceSite& site, int64_t idx);

} // namespace openmc

#endif // OPENMC_PHYSICS_H