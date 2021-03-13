#pragma once

#include "openmc/event.h"
#include "openmc/particle.h"

#include <cub/block/block_scan.cuh>

namespace openmc {
namespace gpu {

/*
 * This does exactly what the following few lines of code describe, but in a
 * more efficient way that avoids atomic operation contentions. The idea is to
 * construct the indices over the current thread block we need to write into
 * with a parallel scan over boolean values regarding whether the fuel or
 * nonfuel queue should be appended to.
 *
 * (TODO put back atomicInc code)
 *
 * This function conditionally appends to one of two arrays living on the GPU
 * based on the values of into_q1 and into_q2. There are two queues in this
 * function, denoted by the abbreviation "q1" and "q2".  The values of into_q1
 * and into_q2 denote whether the particle provided (represented by the pair
 * (p, p_idx)) should go into either queue 1 or queue 2. Consequently, these
 * booleans _should_ be mutually exclusive, but there's nothing stopping both
 * of them from being true. It is both possible and frequently occurring that
 * into_q1 and into_q2 are both false, in which case the particle under
 * consideration has just died due to absorption or leakage.
 *
 * The use of this function cuts the execution time of the GPU kernels by
 * around 50% compared to the use of atomicInc instructions for writing into
 * the particle queues. The key thing here is that we only need one atomicAdd
 * per thread block rather than one atomicInc per thread.
 */

template<unsigned BLOCK_SIZE>
__device__ void block_queue_pushback(bool const& into_q1, bool const& into_q2,
  EventQueueItem* const& q1, EventQueueItem* const& q2,
  unsigned* const& q1_index, unsigned* const& q2_index, Particle* const& p,
  unsigned const& p_idx)
{
  // Casts bools to integers of either 0 or 1. We will apply
  // a parallel prefix sum to these to transform them into
  // the indices we need to write to.
  unsigned goes_to_q1 = into_q1;
  unsigned goes_to_q2 = into_q2;

  // Boilerplate for block-level scan and reduction
  typedef cub::BlockScan<unsigned, BLOCK_SIZE> BlockScanT;
  __shared__ typename BlockScanT::TempStorage scan;
  __shared__ unsigned start;

  // First put nonfuel indices back to their queue
  BlockScanT(scan).ExclusiveSum(goes_to_q1, goes_to_q1);
  constexpr unsigned max_thread_idx = BLOCK_SIZE - 1;
  if (threadIdx.x == max_thread_idx)
    start = atomicAdd(q1_index, goes_to_q1 + ((unsigned)into_q1));
  __syncthreads();
  if (into_q1)
    q1[start + goes_to_q1] = {*p, p_idx};

  // Now put fuel indices back to their queue
  BlockScanT(scan).ExclusiveSum(goes_to_q2, goes_to_q2);
  if (threadIdx.x == max_thread_idx)
    start = atomicAdd(q2_index, goes_to_q2 + ((unsigned)into_q2));
  __syncthreads();
  if (into_q2)
    q2[start + goes_to_q2] = {*p, p_idx};
}

} // namespace gpu
} // namespace openmc
