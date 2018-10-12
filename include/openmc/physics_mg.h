//! \file physics_mg.h
//! Methods needed to perform the collision physics for multi-group mode

#ifndef OPENMC_PHYSICS_MG_H
#define OPENMC_PHYSICS_MG_H

#include "openmc/capi.h"
#include "openmc/particle.h"
#include "openmc/nuclide.h"

namespace openmc {

//TODO: Remove energy_bin_avg and material_xs parameters when they reside on
// the C-side this should happen after materials, physics, input, and tallies
// are brought over

//! \brief samples particle behavior after a collision event.
extern "C" void
collision_mg(Particle* p, const double* energy_bin_avg,
     const MaterialMacroXS& material_xs);

//! \brief samples a reaction type.
//!
//! Note that there is special logic when suvival biasing is turned on since
//! fission and disappearance are treated implicitly.
void
sample_reaction(Particle* p, const double* energy_bin_avg,
     const MaterialMacroXS& material_xs);

//! \brief Samples the scattering event
void
scatter(Particle* p, const double* energy_bin_avg);

//! \brief Determines the average total, prompt and delayed neutrons produced
//! from fission and creates the appropriate bank sites.
void
create_fission_sites(Particle* p, Bank* bank_array, int64_t& size_bank,
     const int64_t& bank_array_size, const MaterialMacroXS& material_xs);

//! \brief Handles an absorption event
void
absorption(Particle* p, const MaterialMacroXS& material_xs);

} // namespace openmc
#endif // OPENMC_PHYSICS_MG_H
