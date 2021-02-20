#include "openmc/cuda/calculate_xs.h"
#include "openmc/cuda/cross_surface.h"

namespace openmc {
namespace gpu {

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

    // Particle now needs an XS lookup
    if (p.alive_) {
      if (p.material_ == MATERIAL_VOID ||
          !gpu::materials[p.material_]->fissionable_)
        calculate_nonfuel_xs_queue[atomicInc(
          &managed_calculate_nonfuel_queue_index, 0xFFFFFFF)] = {p, p_idx};
      else
        calculate_fuel_xs_queue[atomicInc(
          &managed_calculate_fuel_queue_index, 0xFFFFFFF)] = {p, p_idx};
    }
  }
}

} // namespace gpu
} // namespace openmc
