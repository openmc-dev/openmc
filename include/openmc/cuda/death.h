#pragma once

#include "openmc/particle.h"
#include "openmc/settings.h" // BLOCKSIZE

namespace openmc {
namespace gpu {

__global__ __launch_bounds__(BLOCKSIZE) void process_death_events_device(
  unsigned n_particles)
{
  unsigned tid = threadIdx.x + blockDim.x * blockIdx.x;
  Particle p(tid);
  if (tid < n_particles)
    p.event_death();
}

} // namespace gpu
} // namespace openmc
