#pragma once

#include "openmc/constants.h" // MATERIAL_VOID
#include "openmc/event.h"
#include "openmc/material.h"
#include "openmc/memory.h"
#include "openmc/nuclide.h"
#include "openmc/particle.h"

namespace openmc {
namespace gpu {

// Constants used for XS lookup often accessed in this kernel, hence
// meriting the use of constant memory
// TODO: put these variables where you'd expect them to be, and
// remove the inclusion of this header where it is no longer necessary.
extern __constant__ unique_ptr<Material>* materials;
extern __constant__ unique_ptr<Nuclide>* nuclides;
extern __constant__ Particle* particles;
extern __constant__ NuclideMicroXS* micros;
extern __constant__ double energy_min_neutron;
extern __constant__ double energy_max_neutron;
extern __constant__ double log_spacing;
extern __constant__ unsigned number_nuclides;
extern __constant__ bool need_depletion_rx;

extern __managed__ unsigned managed_calculate_fuel_queue_index;
extern __managed__ unsigned managed_calculate_nonfuel_queue_index;

__global__ void process_calculate_xs_events_device(
  EventQueueItem* queue, unsigned queue_size);

} // namespace gpu
}
