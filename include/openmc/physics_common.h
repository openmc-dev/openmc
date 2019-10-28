//! \file physics_common.h
//! A collection of physics methods common to MG, CE, photon, etc.

#ifndef OPENMC_PHYSICS_COMMON_H
#define OPENMC_PHYSICS_COMMON_H

#include "openmc/particle.h"

namespace openmc {

//! \brief Performs the russian roulette operation for a particle
extern "C" void
russian_roulette(Particle* p);

} // namespace openmc
#endif // OPENMC_PHYSICS_COMMON_H
