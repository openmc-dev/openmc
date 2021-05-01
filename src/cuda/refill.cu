#include "openmc/cuda/refill.h"

namespace openmc {
namespace gpu {

__managed__ unsigned dead_particle_indices_indx;
__constant__ unsigned* dead_particle_indices;

__global__ void scan_for_dead_particles(unsigned n_particles)
{
  unsigned tid = threadIdx.x + blockDim.x * blockIdx.x;
  Particle p(tid);

  if (tid < n_particles) {
    if (!p.alive()) {
      p.event_death();
      auto loc = atomicInc(&dead_particle_indices_indx, n_particles);
      dead_particle_indices[loc] = tid;
    }
  }
}

} // namespace gpu
} // namespace openmc
