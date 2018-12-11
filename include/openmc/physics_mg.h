//! \file physics_mg.h
//! Methods needed to perform the collision physics for multi-group mode

#ifndef OPENMC_PHYSICS_MG_H
#define OPENMC_PHYSICS_MG_H

#include "openmc/capi.h"
#include "openmc/particle.h"
#include "openmc/nuclide.h"

namespace openmc {

//! \brief samples particle behavior after a collision event.
//! \param p Particle to operate on
extern "C" void
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
//! \param bank_array The particle bank to populate
//! \param size_bank Number of particles currently in the bank
//! \param bank_array_size Allocated size of the bank
void
create_fission_sites(Particle* p, Bank* bank_array, int64_t* size_bank,
     int64_t bank_array_size);

//! \brief Handles an absorption event
//! \param p Particle to operate on
void
absorption(Particle* p);

} // namespace openmc
#endif // OPENMC_PHYSICS_MG_H
