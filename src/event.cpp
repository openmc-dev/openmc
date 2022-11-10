#include "openmc/bank.h"
#include "openmc/cell.h"
#include "openmc/event.h"
#include "openmc/lattice.h"
#include "openmc/material.h"
#include "openmc/simulation.h"
#include "openmc/sort.h"
#include "openmc/surface.h"
#include "openmc/timer.h"
#include "openmc/tallies/tally.h"
#ifndef DEVICE_PRINTF
#define printf(fmt, ...) (0)
#endif

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace simulation {

SharedArray<EventQueueItem> calculate_fuel_xs_queue;
SharedArray<EventQueueItem> calculate_nonfuel_xs_queue;
SharedArray<EventQueueItem> advance_particle_queue;
SharedArray<EventQueueItem> surface_crossing_queue;
SharedArray<EventQueueItem> collision_queue;
SharedArray<EventQueueItem> revival_queue;

int current_source_offset;

int sort_counter{0};

} // namespace simulation

//==============================================================================
// Non-member functions
//==============================================================================

void sort_queue(SharedArray<EventQueueItem>& queue)
{
  simulation::time_event_sort.start();

  if (queue.size() > settings::minimum_sort_items)
  {
    simulation::sort_counter++;

    #ifdef CUDA_THRUST_SORT
    device_sort_event_queue_item(queue.device_data(), queue.device_data() + queue.size());
    #elif SYCL_SORT
    sort_queue_SYCL(queue.device_data(), queue.device_data() + queue.size());
    #else
    // Transfer queue information to the host
    #pragma omp target update from(queue.data_[:queue.size()])

    // Sort queue via OpenMP parallel sort implementation
    quickSort_parallel(queue.data(), queue.size());

    // Transfer queue information back to the device
    #pragma omp target update to(queue.data_[:queue.size()])
    #endif
  }

  simulation::time_event_sort.stop();
}

bool is_sorted(SharedArray<EventQueueItem>& queue)
{
  int not_sorted = 0;
  if( simulation::calculate_fuel_xs_queue.size() > settings::minimum_sort_items )
  {
    #pragma omp target teams distribute parallel for reduction(+:not_sorted)
    for( int i = 1; i < queue.size(); i++ )
    {
      if( queue[i-1].E > queue[i].E )
        not_sorted += 1;
    }
  }
  if( !not_sorted )
    return true;
  else
    return false;
}

void init_event_queues(int n_particles)
{
  simulation::calculate_fuel_xs_queue.reserve(n_particles);
  simulation::calculate_nonfuel_xs_queue.reserve(n_particles);
  simulation::advance_particle_queue.reserve(n_particles);
  simulation::surface_crossing_queue.reserve(n_particles);
  simulation::collision_queue.reserve(n_particles);
  simulation::revival_queue.reserve(n_particles);

  simulation::particles.resize(n_particles);

  // Allocate any queues that are needed on device
  simulation::calculate_fuel_xs_queue.allocate_on_device();
  simulation::calculate_nonfuel_xs_queue.allocate_on_device();
  simulation::advance_particle_queue.allocate_on_device();
  simulation::surface_crossing_queue.allocate_on_device();
  simulation::collision_queue.allocate_on_device();
  simulation::revival_queue.allocate_on_device();

  #pragma omp target update to(simulation::work_per_rank)
}

void free_event_queues(void)
{
  simulation::calculate_fuel_xs_queue.clear();
  simulation::calculate_nonfuel_xs_queue.clear();
  simulation::advance_particle_queue.clear();
  simulation::surface_crossing_queue.clear();
  simulation::collision_queue.clear();
  simulation::revival_queue.clear();

  simulation::particles.clear();
}

void dispatch_xs_event(int buffer_idx)
{
  Particle& p = simulation::device_particles[buffer_idx];

  // Determine if the particle requires an XS lookup or not
  bool needs_lookup = p.event_calculate_xs_dispatch();

  if (needs_lookup) {
    // If a lookup is needed, dispatch to fuel vs. non-fuel lookup queue
    if (!model::materials[p.material_].fissionable_) {
      simulation::calculate_nonfuel_xs_queue.thread_safe_append({p.E_, p.material_, buffer_idx});
    } else {
      simulation::calculate_fuel_xs_queue.thread_safe_append({p.E_, buffer_idx});
    }
  } else {
    // Otherwise, particle can move directly to the advance particle queue
    simulation::advance_particle_queue.thread_safe_append({p.E_, buffer_idx});
  }
}

void process_init_events(int n_particles)
{
  simulation::time_event_init.start();

  simulation::current_source_offset = n_particles;
  #pragma omp target update to(simulation::current_source_offset)

  double total_weight = 0.0;

  #pragma omp target teams distribute parallel for
  for (int i = 0; i < n_particles; i++) {
    initialize_history(simulation::device_particles[i], i + 1);
    dispatch_xs_event(i);
  }

  // The loop below can in theory be combined with the one above,
  // but is present here as a compiler bug workaround
  #pragma omp target teams distribute parallel for reduction(+:total_weight)
  for (int i = 0; i < n_particles; i++) {
    total_weight += simulation::device_particles[i].wgt_;
  }
  simulation::time_event_init.stop();

  // Write total weight to global variable
  simulation::total_weight = total_weight;

  simulation::calculate_fuel_xs_queue.sync_size_device_to_host();
  simulation::calculate_nonfuel_xs_queue.sync_size_device_to_host();
  simulation::advance_particle_queue.sync_size_device_to_host();

}

bool depletion_rx_check()
{
  return !model::active_tracklength_tallies.empty() &&
    (simulation::need_depletion_rx || simulation::depletion_scores_present);
}

void process_calculate_xs_events_nonfuel()
{
  // Sort non fuel lookup queue by material and energy
  sort_queue(simulation::calculate_nonfuel_xs_queue);

  simulation::time_event_calculate_xs.start();
  simulation::time_event_calculate_xs_nonfuel.start();

  bool need_depletion_rx = depletion_rx_check();

  int offset = simulation::advance_particle_queue.size();;

  #pragma omp target teams distribute parallel for
  for (int i = 0; i < simulation::calculate_nonfuel_xs_queue.size(); i++) {
    int buffer_idx = simulation::calculate_nonfuel_xs_queue[i].idx;
    Particle& p = simulation::device_particles[buffer_idx];
    p.event_calculate_xs_execute(need_depletion_rx);
    simulation::advance_particle_queue[offset + i] = simulation::calculate_nonfuel_xs_queue[i];
  }

  // After executing a calculate_xs event, particles will
  // always require an advance event. Therefore, we don't need to use
  // the protected enqueuing function.
  simulation::advance_particle_queue.resize(offset + simulation::calculate_nonfuel_xs_queue.size());
  simulation::calculate_nonfuel_xs_queue.resize(0);

  simulation::time_event_calculate_xs.stop();
  simulation::time_event_calculate_xs_nonfuel.stop();
}

void process_calculate_xs_events_fuel()
{
  // Sort fuel lookup queue by energy
  sort_queue(simulation::calculate_fuel_xs_queue);

  // The below line can be used to check if the queue has actually been sorted.
  // May be useful for debugging future on-device sorting strategies.
  //assert(is_sorted(simulation::calculate_fuel_xs_queue));

  simulation::time_event_calculate_xs.start();
  simulation::time_event_calculate_xs_fuel.start();

  bool need_depletion_rx = depletion_rx_check();

  int offset = simulation::advance_particle_queue.size();;

  #pragma omp target teams distribute parallel for
  for (int i = 0; i < simulation::calculate_fuel_xs_queue.size(); i++) {
    int buffer_idx = simulation::calculate_fuel_xs_queue[i].idx;
    Particle& p = simulation::device_particles[buffer_idx];
    p.event_calculate_xs_execute(need_depletion_rx);
    simulation::advance_particle_queue[offset + i] = simulation::calculate_fuel_xs_queue[i];
  }

  // After executing a calculate_xs event, particles will
  // always require an advance event. Therefore, we don't need to use
  // the protected enqueuing function.
  simulation::advance_particle_queue.resize(offset + simulation::calculate_fuel_xs_queue.size());
  simulation::calculate_fuel_xs_queue.resize(0);

  simulation::time_event_calculate_xs.stop();
  simulation::time_event_calculate_xs_fuel.stop();
}

void process_advance_particle_events(int n_particles)
{
  simulation::time_event_advance_particle.start();

  #pragma omp target teams distribute parallel for
  for (int i = 0; i < simulation::advance_particle_queue.size(); i++) {
    int buffer_idx = simulation::advance_particle_queue[i].idx;
    Particle& p = simulation::device_particles[buffer_idx];
    p.event_advance();
    if (p.collision_distance_ > p.boundary_.distance) {
      simulation::surface_crossing_queue.thread_safe_append({p.E_, buffer_idx});
    } else {
      simulation::collision_queue.thread_safe_append({p.E_, buffer_idx});
    }
  }
  simulation::surface_crossing_queue.sync_size_device_to_host();
  simulation::collision_queue.sync_size_device_to_host();
  
  simulation::time_event_advance_particle.stop();
  simulation::time_event_tally.start();
  
  // Perform tracklength tallying on device if active tallies are present
  if (!model::active_tracklength_tallies.empty()) {
    bool need_depletion_rx = depletion_rx_check();
    #pragma omp target teams distribute parallel for
    for (int i = 0; i < simulation::advance_particle_queue.size(); i++) {
      int buffer_idx = simulation::advance_particle_queue[i].idx;
      Particle& p = simulation::device_particles[buffer_idx];
      p.event_tracklength_tally(need_depletion_rx);
    }
  }

  simulation::advance_particle_queue.resize(0);
  simulation::time_event_tally.stop();
}

void process_surface_crossing_events()
{
  simulation::time_event_surface_crossing.start();

  #pragma omp target teams distribute parallel for
  for (int i = 0; i < simulation::surface_crossing_queue.size(); i++) {
    int buffer_idx = simulation::surface_crossing_queue[i].idx;
    Particle& p = simulation::device_particles[buffer_idx];
    p.event_cross_surface();
    if (p.alive())
      dispatch_xs_event(buffer_idx);
    else
      simulation::revival_queue.thread_safe_append({p.E_, buffer_idx});
  }

  simulation::calculate_fuel_xs_queue.sync_size_device_to_host();
  simulation::calculate_nonfuel_xs_queue.sync_size_device_to_host();
  simulation::advance_particle_queue.sync_size_device_to_host();
  simulation::revival_queue.sync_size_device_to_host();
  simulation::surface_crossing_queue.resize(0);

  simulation::time_event_surface_crossing.stop();
}

void process_collision_events()
{
  simulation::time_event_collision.start();

  #pragma omp target teams distribute parallel for
  for (int i = 0; i < simulation::collision_queue.size(); i++) {
    int buffer_idx = simulation::collision_queue[i].idx;
    Particle& p = simulation::device_particles[buffer_idx];
    p.event_collide();
    if (p.alive())
      dispatch_xs_event(buffer_idx);
    else
      simulation::revival_queue.thread_safe_append({p.E_, buffer_idx});
  }

  simulation::calculate_fuel_xs_queue.sync_size_device_to_host();
  simulation::calculate_nonfuel_xs_queue.sync_size_device_to_host();
  simulation::advance_particle_queue.sync_size_device_to_host();
  simulation::revival_queue.sync_size_device_to_host();
  simulation::collision_queue.resize(0);

  simulation::time_event_collision.stop();
}


void process_death_events(int n_particles)
{
  simulation::time_event_death.start();

  // Local keff tally accumulators
  // (workaround for compilers that don't like reductions w/global variables)
  double absorption = 0.0;
  double collision = 0.0;
  double tracklength = 0.0;
  double leakage = 0.0;

  #pragma omp target teams distribute parallel for reduction(+:absorption, collision, tracklength, leakage)
  for (int i = 0; i < n_particles; i++) {
    Particle& p = simulation::device_particles[i];
    p.accumulate_keff_tallies_local(absorption, collision, tracklength, leakage);
  }

  // Write local reduction results to global values
  global_tally_absorption  += absorption;
  global_tally_collision   += collision;
  global_tally_tracklength += tracklength;
  global_tally_leakage     += leakage;

  // Move particle progeny count array back to host
  #pragma omp target update from(simulation::device_progeny_per_particle[:simulation::progeny_per_particle.size()])

  simulation::time_event_death.stop();

}

void process_revival_events()
{
  simulation::time_event_revival.start();

  // Accumulator for particle weights from any sourced particles
  double extra_weight = 0;

  #pragma omp target teams distribute parallel for reduction(+:extra_weight)
  for (int i = 0; i < simulation::revival_queue.size(); i++) {
    int buffer_idx = simulation::revival_queue[i].idx;
    Particle& p = simulation::device_particles[buffer_idx];

    // First, attempt to revive a secondary particle
    p.event_revive_from_secondary();

    // If no secondary particle was found, attempt to source new particle
    if (!p.alive()) {

      // Before sourcing new particle, execute death event for old particle
      p.event_death();

      // Atomically retrieve source offset for new particle
      int source_offset_idx;
      #pragma omp atomic capture //seq_cst
      source_offset_idx = simulation::current_source_offset++;

      // Check that we are not going to run more particles than the user specified
      if (source_offset_idx < simulation::work_per_rank) {
        // If a valid particle is sourced, initialize it and accumulate its weight
        initialize_history(p, source_offset_idx + 1);
        extra_weight += p.wgt_;
      }
    }
    // If particle has either been revived or a new particle has been sourced,
    // then dispatch particle to appropriate queue. Otherwise, if particle is
    // dead, then it will not be queued anywhere.
    if (p.alive())
      dispatch_xs_event(buffer_idx);
  }

  simulation::revival_queue.resize(0);
  simulation::calculate_fuel_xs_queue.sync_size_device_to_host();
  simulation::calculate_nonfuel_xs_queue.sync_size_device_to_host();
  simulation::advance_particle_queue.sync_size_device_to_host();

  // Add any newly sourced particle weights to global variable
  simulation::total_weight += extra_weight;

  simulation::time_event_revival.stop();
}

} // namespace openmc
