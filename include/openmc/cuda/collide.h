#pragma once

#include "openmc/event.h"
#include "openmc/particle.h"

namespace openmc {
namespace gpu {

__global__ void process_collision_events_device(EventQueueItem* queue,
  unsigned queue_size, EventQueueItem* calculate_nonfuel_xs_queue,
  EventQueueItem* calculate_fuel_xs_queue);

} // namespace gpu
} // namespace openmc
