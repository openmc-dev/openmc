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
namespace gpu {

// Constants used for XS lookup often accessed in this kernel, hence
// meriting the use of constant memory
extern __constant__ unique_ptr<Material>* materials;
extern __constant__ unique_ptr<Nuclide>* nuclides;
extern __constant__ Particle* particles;
extern __constant__ double energy_min_neutron;
extern __constant__ double log_spacing;
extern __constant__ bool need_depletion_rx;

__global__ void process_calculate_xs_events_device(
  EventQueueItem* queue, unsigned queue_size);

} // namespace gpu
}
