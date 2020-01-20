//! \file physics_mg.h
//! Methods needed to perform the collision physics for multi-group mode

#ifndef OPENMC_PHYSICS_MG_H
#define OPENMC_PHYSICS_MG_H

#include "openmc/capi.h"
#include "openmc/particle.h"
#include "openmc/nuclide.h"

#include <vector>

namespace openmc {

//! \brief samples particle behavior after a collision event.
//! \param p Particle to operate on
void
collision_mg(Particle* p);

//! \brief samples a reaction type.
//!
//! Note that there is special logic when suvival biasing is turned on since
//! fission and disappearance are treated implicitly.
//! \param p Particle to operate on
void
sample_reaction(Particle* p);

//! \brief Samples the scattering event
//! \param p Particle to operate on
void
scatter(Particle* p);

//! \brief Determines the average total, prompt and delayed neutrons produced
//! from fission and creates the appropriate bank sites.
//! \param p Particle to operate on
void
create_fission_sites(Particle* p);

//! \brief Handles an absorption event
//! \param p Particle to operate on
void
absorption(Particle* p);

} // namespace openmc
#endif // OPENMC_PHYSICS_MG_H
