//! \file physics_mg.h
//! Methods needed to perform the collision physics for multi-group mode

#ifndef OPENMC_PHYSICS_MG_H
#define OPENMC_PHYSICS_MG_H

#include "openmc/particle.h"

namespace openmc {

//==============================================================================
// SCATTER
//==============================================================================

//! \brief Samples the scattering event
extern "C" void
scatter(Particle* p, const double* energy_bin_avg);

} // namespace openmc
#endif // OPENMC_PHYSICS_MG_H
