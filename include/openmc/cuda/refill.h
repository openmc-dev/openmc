#pragma once

#include "openmc/cuda/block_queue_pushback.h"
#include "openmc/cuda/calculate_xs.h"
#include "openmc/particle.h"
#include "openmc/simulation.h" // initialize_history

namespace openmc {
namespace gpu {

extern __managed__ unsigned dead_particle_indices_indx;
extern __constant__ unsigned* dead_particle_indices;

__global__ void scan_for_dead_particles(unsigned n_particles);

template<unsigned BLOCK_SIZE>
__global__ void refill_dead_particle_slots(unsigned n_refilled,
  unsigned source_offset,
  EventQueueItem* __restrict__ calculate_nonfuel_xs_queue,
  EventQueueItem* __restrict__ calculate_fuel_xs_queue)
{
  unsigned tid = threadIdx.x + blockDim.x * blockIdx.x;
  bool nonfuel = false;
  bool fuel = false;
  unsigned p_idx = tid < n_refilled ? dead_particle_indices[tid] : 0;
  Particle p(p_idx);

  if (tid < n_refilled) {
    p.initialize_values();
    initialize_history(p, source_offset + tid + 1);
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
