#include "openmc/cuda/advance.h"
#include "openmc/cuda/calculate_xs.h" // for particle array pointer in constant memory
#include "openmc/geometry.h"
#include "openmc/random_lcg.h"

namespace openmc {
namespace gpu {
__managed__ unsigned managed_surface_crossing_queue_index;
__managed__ unsigned managed_collision_queue_index;

__global__ void process_advance_events_device(EventQueueItem* queue,
  unsigned queue_size, EventQueueItem* surface_crossing_queue,
  EventQueueItem* collision_queue)
{
  for (unsigned tid = threadIdx.x + blockDim.x * blockIdx.x; tid < queue_size;
       tid += blockDim.x * gridDim.x) {
    Particle& p = particles[queue[tid].idx];
    p.boundary_ = distance_to_boundary(p);

    // Assuming only tracking neutrons on the GPU for now
    double distance = -log(prn(p.seeds_)) / p.macro_xs_.total;
    distance = min(distance, p.boundary_.distance);

    // Advance the particle
    for (int j = 0; j < p.n_coord_; ++j)
      p.coord_[j].r += distance * p.coord_[j].u;

    // TODO add tracklength tallies here

    // Always score the the k_eff tracklength tally. It's simply not
    // used if doing a fixed source calculation.
    p.keff_tally_tracklength_ += p.wgt_ * distance * p.macro_xs_.nu_fission;

    // TODO differential track tally stuff

    // Push back to main arrays
    if (distance == p.boundary_.distance)
      surface_crossing_queue[atomicInc(&managed_surface_crossing_queue_index,
        0xFFFFFFF)] = {p, queue[tid].idx};
    else
      collision_queue[atomicInc(&managed_collision_queue_index, 0xFFFFFFF)] = {
        p, queue[tid].idx};
  }
}

} // namespace gpu
} // namespace openmc
