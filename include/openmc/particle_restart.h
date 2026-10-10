#ifndef OPENMC_PARTICLE_RESTART_H
#define OPENMC_PARTICLE_RESTART_H

#include "openmc/particle.h"

namespace openmc {

void run_particle_restart();
void run_lost_particle_track(Particle& lost);

} // namespace openmc

#endif // OPENMC_PARTICLE_RESTART_H
