#pragma once

#include "openmc/event.h"
#include "openmc/particle.h"

#include "openmc/cuda/block_queue_pushback.h"
#include "openmc/cuda/calculate_xs.h"
#include "openmc/cuda/collide.h"

namespace openmc {
namespace gpu {

template<unsigned BLOCK_SIZE>
__global__ void process_collision_events_device(EventQueueItem* queue,
  unsigned queue_size, EventQueueItem* calculate_nonfuel_xs_queue,
  EventQueueItem* calculate_fuel_xs_queue)
{
  unsigned tid = threadIdx.x + blockDim.x * blockIdx.x;
  bool nonfuel = false;
  bool fuel = false;
  unsigned p_idx;
  Particle* p;

  if (tid < queue_size) {
    p_idx = queue[tid].idx;
    p = particles + p_idx;
    p->event_collide();

    // Replace with revival from secondaries eventually
    p->n_event_++;

    // These are used as booleans here, but are converted to indices shortly.
    nonfuel = p->alive_ && (p->material_ == MATERIAL_VOID ||
                             !gpu::materials[p->material_]->fissionable_);
    fuel = p->alive_ && !nonfuel;
  }

  // Particle now needs an XS lookup
  block_queue_pushback<BLOCK_SIZE>(nonfuel, fuel, calculate_nonfuel_xs_queue,
    calculate_fuel_xs_queue, &managed_calculate_nonfuel_queue_index,
    &managed_calculate_fuel_queue_index, p, p_idx);
}

} // namespace gpu
} // namespace openmc
