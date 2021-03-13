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
__global__ void process_advance_events_device(EventQueueItem* queue,
  unsigned queue_size, EventQueueItem* surface_crossing_queue,
  EventQueueItem* collision_queue)
{
  unsigned tid = threadIdx.x + blockDim.x * blockIdx.x;
  bool surface = false;
  bool collision = false;
  unsigned p_idx;
  Particle* p;

  if (tid < queue_size) {
    p_idx = queue[tid].idx;
    p = particles + p_idx;
    p->event_advance();

    surface = p->collision_distance_ > p->boundary_.distance;
    collision = !surface;
  }

  // Push back to either surface crossing or event queue arrays
  block_queue_pushback<BLOCK_SIZE>(surface, collision, surface_crossing_queue,
    collision_queue, &managed_surface_crossing_queue_index,
    &managed_collision_queue_index, p, p_idx);
}

} // namespace gpu
} // namespace openmc
