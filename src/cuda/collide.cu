#include "openmc/cuda/calculate_xs.h"
#include "openmc/cuda/collide.h"

namespace openmc {
namespace gpu {

__global__ void process_collision_events_device(EventQueueItem* queue,
  unsigned queue_size, EventQueueItem* calculate_nonfuel_xs_queue,
  EventQueueItem* calculate_fuel_xs_queue)
{
  for (unsigned tid = threadIdx.x + blockDim.x * blockIdx.x; tid < queue_size;
       tid += blockDim.x * gridDim.x) {
    Particle& p = particles[queue[tid].idx];
    p.event_collide();
    if (p.alive_)
      dispatch_xs_event_device(
        p, queue[tid].idx, calculate_nonfuel_xs_queue, calculate_fuel_xs_queue);
  }
}

} // namespace gpu
} // namespace openmc
