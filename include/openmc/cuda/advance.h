#pragma once

#include "openmc/cuda/block_queue_pushback.h"
#include "openmc/cuda/calculate_xs.h" // gpu::particles
#include "openmc/event.h"
#include "openmc/geometry.h"
#include "openmc/particle.h"
#include "openmc/random_lcg.h"

namespace openmc {
namespace gpu {

extern __managed__ unsigned managed_surface_crossing_queue_index;
extern __managed__ unsigned managed_collision_queue_index;

template<unsigned BLOCK_SIZE>
__global__ void __launch_bounds__(BLOCK_SIZE)
  process_advance_events_device(EventQueueItem* __restrict__ queue,
    unsigned queue_size, EventQueueItem* __restrict__ surface_crossing_queue,
    EventQueueItem* __restrict__ collision_queue)
{
  unsigned tid = threadIdx.x + blockDim.x * blockIdx.x;
  bool surface = false;
  bool collision = false;
  unsigned p_idx = tid < queue_size ? queue[tid].idx : 0;
  Particle p(p_idx);

  if (tid < queue_size) {
    p.event_advance();

    surface = p.collision_distance() > p.boundary().distance;
    collision = !surface;
  }

  // Push back to either surface crossing or event queue arrays
  block_queue_pushback<BLOCK_SIZE>(surface, collision, surface_crossing_queue,
    collision_queue, &managed_surface_crossing_queue_index,
    &managed_collision_queue_index, p, p_idx);
}

} // namespace gpu
} // namespace openmc
