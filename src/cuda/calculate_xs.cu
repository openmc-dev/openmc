#include "openmc/cuda/calculate_xs.h"
namespace openmc {

__global__ void process_calculate_xs_events_device(Particle* particles, EventQueueItem* queue, unsigned queue_size, unique_ptr<Material>* materials, unique_ptr<Nuclide>* nuclides) {
  for (unsigned tid=threadIdx.x+blockDim.x*blockIdx.x; tid<queue_size; tid += blockDim.x * gridDim.x) {
    Particle& p = particles[queue[tid].idx];
    // p.event_calculate_xs();
  }
}

}
