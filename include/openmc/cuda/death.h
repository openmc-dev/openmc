#pragma once

#include "openmc/particle.h"

namespace openmc {
namespace gpu {

__global__ void process_death_events_device(unsigned n_particles)
{
  unsigned tid = threadIdx.x + blockDim.x * blockIdx.x;
  Particle p(tid);
  if (tid < n_particles)
    p.event_death();
}

} // namespace gpu
} // namespace openmc
