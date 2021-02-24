#pragma once

#include "openmc/event.h"
#include "openmc/particle.h"

#include "openmc/cuda/block_queue_pushback.h"
#include "openmc/cuda/calculate_xs.h"
#include "openmc/cuda/cross_surface.h"

namespace openmc {
namespace gpu {

template<unsigned BLOCK_SIZE>
__global__ void process_surface_crossing_events_device(EventQueueItem* queue,
  unsigned queue_size, EventQueueItem* calculate_nonfuel_xs_queue,
  EventQueueItem* calculate_fuel_xs_queue)
{
  unsigned tid = threadIdx.x + blockDim.x * blockIdx.x;
  unsigned goes_to_nonfuel_queue = 0;
  unsigned goes_to_fuel_queue = 0;
  bool fuel = false;
  bool nonfuel = false;
  int64_t p_idx;
  Particle* p;

  if (tid < queue_size) {
    p_idx = queue[tid].idx;
    p = particles + p_idx;
    p->event_cross_surface();

    // Replace with revival from secondaries eventually
    p->n_event_++;

    // These are used as booleans here, but are converted to indices shortly.
    goes_to_nonfuel_queue =
      p->alive_ && (p->material_ == MATERIAL_VOID ||
                     !gpu::materials[p->material_]->fissionable_);
    goes_to_fuel_queue = p->alive_ && !goes_to_nonfuel_queue;
    fuel = goes_to_fuel_queue;
    nonfuel = goes_to_nonfuel_queue;
  }

  block_queue_pushback<BLOCK_SIZE>(nonfuel, fuel, calculate_nonfuel_xs_queue,
    calculate_fuel_xs_queue, &managed_calculate_nonfuel_queue_index,
    &managed_calculate_fuel_queue_index, p, p_idx);
}

} // namespace gpu
} // namespace openmc
