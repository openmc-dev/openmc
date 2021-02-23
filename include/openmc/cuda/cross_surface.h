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
  for (unsigned tid = threadIdx.x + blockDim.x * blockIdx.x; tid < queue_size;
       tid += blockDim.x * gridDim.x) {
    auto const& p_idx = queue[tid].idx;
    Particle& p = particles[p_idx];
    p.event_cross_surface();

    // Replace with revival from secondaries eventually
    p.n_event_++;

    /*
     * The code below here does exactly what the following few lines of code
     * describe, but in a more efficient way that avoids atomic operation
     * contentions. The idea is to construct the indices over the current
     * thread block we need to write into with a parallel scan over boolean
     * values regarding whether the fuel or nonfuel queue should be appended
     * to.
     *
     * if (p.alive_) {
     *   if (p.material_ == MATERIAL_VOID ||
     *       !gpu::materials[p.material_]->fissionable_)
     *     calculate_nonfuel_xs_queue[atomicInc(
     *       &managed_calculate_nonfuel_queue_index, 0xFFFFFFF)] = {p, p_idx};
     *   else
     *     calculate_fuel_xs_queue[atomicInc(
     *       &managed_calculate_fuel_queue_index, 0xFFFFFFF)] = {p, p_idx};
     * }
     */

    // Particle now needs an XS lookup. Below is a lock-free (at the block
    // level) method for writing into the XS lookup queues.
    bool goes_to_nonfuel_queue = p.alive_;
    bool goes_to_fuel_queue = goes_to_nonfuel_queue;
    goes_to_nonfuel_queue =
      goes_to_nonfuel_queue && (p.material_ == MATERIAL_VOID ||
                                 !gpu::materials[p.material_]->fissionable_);
    goes_to_fuel_queue = goes_to_fuel_queue && !goes_to_nonfuel_queue;

    // Boilerplate for block-level scan and reduction
    typedef cub::BlockScan<unsigned, BLOCK_SIZE, cub::BLOCK_SCAN_RAKING>
      BlockScanT;
    typename BlockScanT::TempStorage scan_temp_storage;

    __shared__ unsigned start;
    unsigned my_index;      // position in queue to write to
    unsigned amount_to_add; // result of the summation

    // First put nonfuel indices back to their queue
    BlockScanT(scan_temp_storage)
      .ExclusiveSum(goes_to_nonfuel_queue, my_index, amount_to_add);
    if (threadIdx.x == 0)
      start = atomicAdd(&managed_calculate_nonfuel_queue_index, amount_to_add);
    if (goes_to_nonfuel_queue)
      calculate_nonfuel_xs_queue[start + my_index] = {p, p_idx};

    // Now put fuel indices back to their queue
    BlockScanT(scan_temp_storage)
      .ExclusiveSum(goes_to_fuel_queue, my_index, amount_to_add);
    if (threadIdx.x == 0)
      start = atomicAdd(&managed_calculate_fuel_queue_index, amount_to_add);
    if (goes_to_fuel_queue)
      calculate_fuel_xs_queue[start + my_index] = {p, p_idx};
  }
}

} // namespace gpu
} // namespace openmc
