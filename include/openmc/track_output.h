#ifndef OPENMC_TRACK_OUTPUT_H
#define OPENMC_TRACK_OUTPUT_H

#include "openmc/particle.h"

namespace openmc {

//==============================================================================
// Non-member functions
//==============================================================================

void add_particle_track(Particle& p);
void write_particle_track(Particle& p);
void finalize_particle_track(Particle& p);

} // namespace openmc

#endif // OPENMC_TRACK_OUTPUT_H
