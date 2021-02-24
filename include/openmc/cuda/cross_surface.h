#pragma once

#include "openmc/event.h"
#include "openmc/particle.h"

#include "openmc/cuda/calculate_xs.h"
#include "openmc/cuda/cross_surface.h"
#include <cub/block/block_scan.cuh>

namespace openmc {
namespace gpu {

template<unsigned BLOCK_SIZE>
__global__ void process_surface_crossing_events_device(EventQueueItem* queue,
  unsigned queue_size, EventQueueItem* calculate_nonfuel_xs_queue,
  EventQueueItem* calculate_fuel_xs_queue)
{
  unsigned tid = threadIdx.x + blockDim.x * blockIdx.x;
  bool thread_live = tid < queue_size;
  unsigned goes_to_nonfuel_queue = 0;
  unsigned goes_to_fuel_queue = 0;
  bool fuel = false;
  bool nonfuel = false;
  decltype(queue[0].idx) p_idx;
  Particle* p;

  if (thread_live) {
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

  // Boilerplate for block-level scan and reduction
  typedef cub::BlockScan<unsigned, BLOCK_SIZE> BlockScanT;
  __shared__ typename BlockScanT::TempStorage scan;
  __shared__ unsigned start;

  // First put nonfuel indices back to their queue
  BlockScanT(scan).ExclusiveSum(goes_to_nonfuel_queue, goes_to_nonfuel_queue);
  constexpr unsigned max_thread_idx = BLOCK_SIZE - 1;
  if (threadIdx.x == max_thread_idx)
    start = atomicAdd(&managed_calculate_nonfuel_queue_index,
      goes_to_nonfuel_queue + ((unsigned)nonfuel));
  __syncthreads();
  if (nonfuel)
    calculate_nonfuel_xs_queue[start + goes_to_nonfuel_queue] = {*p, p_idx};

  // Now put fuel indices back to their queue
  BlockScanT(scan).ExclusiveSum(goes_to_fuel_queue, goes_to_fuel_queue);
  if (threadIdx.x == max_thread_idx)
    start = atomicAdd(&managed_calculate_fuel_queue_index,
      goes_to_fuel_queue + ((unsigned)fuel));
  __syncthreads();
  if (fuel)
    calculate_fuel_xs_queue[start + goes_to_fuel_queue] = {*p, p_idx};
}

} // namespace gpu
} // namespace openmc
