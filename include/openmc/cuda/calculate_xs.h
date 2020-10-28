#pragma once

#include "openmc/event.h"
#include "openmc/particle.h"
#include "openmc/material.h"
#include "openmc/memory.h"
#include "openmc/nuclide.h"

//==============================================================================
// CUDA kernels for events
//==============================================================================
namespace openmc {
__global__ void process_calculate_xs_events_device(Particle* particles, EventQueueItem* queue, unsigned queue_size, unique_ptr<Material>* materials, unique_ptr<Nuclide>* nuclides);
}
