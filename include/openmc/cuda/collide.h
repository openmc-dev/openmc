#pragma once

#include "openmc/event.h"
#include "openmc/particle.h"

#include "openmc/cuda/block_queue_pushback.h"
#include "openmc/cuda/calculate_xs.h"

namespace openmc {
namespace gpu {

template<unsigned BLOCK_SIZE>
__global__ void __launch_bounds__(BLOCK_SIZE) process_collision_events_device(
  EventQueueItem* __restrict__ queue, unsigned queue_size,
  EventQueueItem* __restrict__ calculate_nonfuel_xs_queue,
  EventQueueItem* __restrict__ calculate_fuel_xs_queue)
{
  unsigned tid = threadIdx.x + blockDim.x * blockIdx.x;
  bool nonfuel = false;
  bool fuel = false;
  unsigned p_idx = tid < queue_size ? queue[tid].idx : 0;
  Particle p(p_idx);

  if (tid < queue_size) {
    p.event_collide();

    // Replace with revival from secondaries eventually
    p.n_event()++;

    // These are used as booleans here, but are converted to indices shortly.
    nonfuel = p.alive() && (p.material() == MATERIAL_VOID ||
                             !gpu::materials[p.material()]->fissionable_);
    fuel = p.alive() && !nonfuel;
  }

  // Particle now needs an XS lookup
  block_queue_pushback<BLOCK_SIZE>(nonfuel, fuel, calculate_nonfuel_xs_queue,
    calculate_fuel_xs_queue, &managed_calculate_nonfuel_queue_index,
    &managed_calculate_fuel_queue_index, p, p_idx);
}

} // namespace gpu
} // namespace openmc
