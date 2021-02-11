#include "openmc/cuda/advance.h"
#include "openmc/cuda/calculate_xs.h" // for particle array pointer in constant memory
#include "openmc/geometry.h"
#include "openmc/random_lcg.h"

namespace openmc {
namespace gpu {
__managed__ unsigned managed_surface_crossing_queue_index;
__managed__ unsigned managed_collision_queue_index;

__global__ void process_advance_events_device(EventQueueItem* queue,
  unsigned queue_size, EventQueueItem* surface_crossing_queue,
  EventQueueItem* collision_queue)
{
  for (unsigned tid = threadIdx.x + blockDim.x * blockIdx.x; tid < queue_size;
       tid += blockDim.x * gridDim.x) {
    Particle& p = particles[queue[tid].idx];
    p.event_advance();

    // Push back to either surface crossing or event queue arrays
    if (p.collision_distance_ > p.boundary_.distance)
      surface_crossing_queue[atomicInc(&managed_surface_crossing_queue_index,
        0xFFFFFFF)] = {p, queue[tid].idx};
    else
      collision_queue[atomicInc(&managed_collision_queue_index, 0xFFFFFFF)] = {
        p, queue[tid].idx};
  }
}

} // namespace gpu
} // namespace openmc
