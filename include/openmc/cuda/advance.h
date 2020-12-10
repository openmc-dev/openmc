#pragma once

#include "openmc/event.h"
#include "openmc/particle.h"

namespace openmc {
namespace gpu {

extern __managed__ unsigned managed_surface_crossing_queue_index;
extern __managed__ unsigned managed_collision_queue_index;
__global__ void process_advance_events_device(EventQueueItem* queue,
  unsigned queue_size, EventQueueItem* surface_crossing_queue,
  EventQueueItem* collision_queue);

} // namespace gpu
} // namespace openmc
