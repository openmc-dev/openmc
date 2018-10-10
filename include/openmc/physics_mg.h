//! \file physics_mg.h
//! Methods needed to perform the collision physics for multi-group mode

#ifndef OPENMC_PHYSICS_MG_H
#define OPENMC_PHYSICS_MG_H

#include "openmc/capi.h"
#include "openmc/particle.h"
#include "openmc/nuclide.h"

namespace openmc {

//==============================================================================
// SCATTER
//==============================================================================

//! \brief Samples the scattering event
extern "C" void
scatter(Particle* p, const double* energy_bin_avg);

//! \brief Determines the average total, prompt and delayed neutrons produced
//! from fission and creates the appropriate bank sites.
extern "C" void
create_fission_sites(Particle* p, Bank* bank_array, int64_t& size_bank,
                     int64_t& bank_array_size, MaterialMacroXS& material_xs);

//! \brief Handles an absorption event
extern "C" void
absorption(Particle* p, MaterialMacroXS& material_xs);

} // namespace openmc
#endif // OPENMC_PHYSICS_MG_H
