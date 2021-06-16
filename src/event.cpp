#include "openmc/bank.h"
#include "openmc/cell.h"
#include "openmc/event.h"
#include "openmc/lattice.h"
#include "openmc/material.h"
#include "openmc/simulation.h"
#include "openmc/surface.h"
#include "openmc/timer.h"
#include "openmc/tallies/tally.h"

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

std::vector<Particle>  particles;
Particle*  device_particles;

} // namespace simulation

//==============================================================================
// Non-member functions
//==============================================================================

void init_event_queues(int64_t n_particles)
{
  simulation::calculate_fuel_xs_queue.reserve(n_particles);
  simulation::calculate_nonfuel_xs_queue.reserve(n_particles);
  simulation::advance_particle_queue.reserve(n_particles);
  simulation::surface_crossing_queue.reserve(n_particles);
  simulation::collision_queue.reserve(n_particles);

  simulation::particles.resize(n_particles);

  // Allocate any queues that are needed on device
  simulation::calculate_fuel_xs_queue.allocate_on_device();
  simulation::calculate_nonfuel_xs_queue.allocate_on_device();
  simulation::advance_particle_queue.allocate_on_device();
  simulation::surface_crossing_queue.allocate_on_device();
  simulation::collision_queue.allocate_on_device();
}

void free_event_queues(void)
{
  simulation::calculate_fuel_xs_queue.clear();
  simulation::calculate_nonfuel_xs_queue.clear();
  simulation::advance_particle_queue.clear();
  simulation::surface_crossing_queue.clear();
  simulation::collision_queue.clear();

  simulation::particles.clear();
}

void dispatch_xs_event(int64_t buffer_idx)
{
  Particle& p = simulation::device_particles[buffer_idx];
  if (p.material_ == MATERIAL_VOID || !model::materials[p.material_].fissionable_) {
    simulation::calculate_nonfuel_xs_queue.thread_safe_append({p, buffer_idx});
  } else {
    simulation::calculate_fuel_xs_queue.thread_safe_append({p, buffer_idx});
  }
}

void process_init_events(int64_t n_particles, int64_t source_offset)
{
  simulation::time_event_init.start();
  
  double total_weight = 0.0;

  #ifdef USE_DEVICE
  #pragma omp target teams distribute parallel for
  #else
  #pragma omp parallel for schedule(runtime)
  #endif
  for (int64_t i = 0; i < n_particles; i++) {
    initialize_history(simulation::device_particles[i], source_offset + i + 1);
    dispatch_xs_event(i);
  }

  // The loop below can in theory be combined with the one above,
  // but is present here as a compiler bug workaround
  #ifdef USE_DEVICE
  // This pragma results in illegal memory errors at runtime
  //#pragma omp target teams distribute parallel for reduction(+:total_weight)
  // This pragma works but runs in serial on device
  #pragma omp target reduction(+:total_weight)
  #else
  #pragma omp parallel for schedule(runtime) reduction(+:total_weight)
  #endif
  for (int64_t i = 0; i < n_particles; i++) {
    total_weight += simulation::device_particles[i].wgt_;
  }

  // Write total weight to global variable
  simulation::total_weight = total_weight;
  
  simulation::calculate_fuel_xs_queue.sync_size_device_to_host();
  simulation::calculate_nonfuel_xs_queue.sync_size_device_to_host();

  simulation::time_event_init.stop();
}

void process_calculate_xs_events(SharedArray<EventQueueItem>& queue)
{
  simulation::time_event_calculate_xs.start();

  // TODO: If using C++17, perform a parallel sort of the queue
  // by particle type, material type, and then energy, in order to
  // improve cache locality and reduce thread divergence on GPU. Prior
  // to C++17, std::sort is a serial only operation, which in this case
  // makes it too slow to be practical for most test problems.
  //
  // std::sort(std::execution::par_unseq, queue.data(), queue.data() + queue.size());

  int64_t offset = simulation::advance_particle_queue.size();;
  
  #ifdef USE_DEVICE
  // This pragma results in illegal memory errors at runtime
  //#pragma omp target teams distribute parallel for
  // This pragma works but only runs one thread per block
  #pragma omp target teams distribute
  #else
  #pragma omp parallel for schedule(runtime)
  #endif
  for (int64_t i = 0; i < queue.size(); i++) {
    int64_t buffer_idx = queue[i].idx;
    Particle& p = simulation::device_particles[buffer_idx];
    p.event_calculate_xs();
    simulation::advance_particle_queue[offset + i] = queue[i];
  }
  
  // After executing a calculate_xs event, particles will
  // always require an advance event. Therefore, we don't need to use
  // the protected enqueuing function.
  simulation::advance_particle_queue.resize(offset + queue.size());

  queue.resize(0);

  simulation::time_event_calculate_xs.stop();
}

void process_advance_particle_events()
{
  simulation::time_event_advance_particle.start();

  #ifdef USE_DEVICE
  #pragma omp target teams distribute parallel for
  #else
  #pragma omp parallel for schedule(runtime)
  #endif
  for (int64_t i = 0; i < simulation::advance_particle_queue.size(); i++) {
    int64_t buffer_idx = simulation::advance_particle_queue[i].idx;
    Particle& p = simulation::device_particles[buffer_idx];
    p.event_advance();
    p.event_advance_tally();

    if (p.collision_distance_ > p.boundary_.distance) {
      simulation::surface_crossing_queue.thread_safe_append({p, buffer_idx});
    } else {
      simulation::collision_queue.thread_safe_append({p, buffer_idx});
    }
  }
  simulation::surface_crossing_queue.sync_size_device_to_host();
  simulation::collision_queue.sync_size_device_to_host();
  simulation::advance_particle_queue.resize(0);

  // DEBUGGING REGION
  /*
  // Let's look at particle 1
  for (int64_t i = 0; i < simulation::advance_particle_queue.size(); i++) {
    int64_t buffer_idx = simulation::advance_particle_queue[i].idx;
    Particle& p = simulation::particles[buffer_idx];
    if( p.id_ == 1 )
    {
      printf("host   particle %ld distance: %.4le\n", p.id_, p.advance_distance_);
      printf("host   particle %ld material: %d\n", p.id_, p.material_);
    }
  }
  exit(1);
  */

  simulation::time_event_advance_particle.stop();
}

void process_surface_crossing_events()
{
  simulation::time_event_surface_crossing.start();

  #ifdef USE_DEVICE
  #pragma omp target teams distribute parallel for
  #else
  #pragma omp parallel for schedule(runtime)
  #endif
  for (int64_t i = 0; i < simulation::surface_crossing_queue.size(); i++) {
    int64_t buffer_idx = simulation::surface_crossing_queue[i].idx;
    Particle& p = simulation::device_particles[buffer_idx];
    p.event_cross_surface();
    p.event_revive_from_secondary();
    if (p.alive_)
      dispatch_xs_event(buffer_idx);
  }

  simulation::calculate_fuel_xs_queue.sync_size_device_to_host();
  simulation::calculate_nonfuel_xs_queue.sync_size_device_to_host();
  simulation::surface_crossing_queue.resize(0);

  simulation::time_event_surface_crossing.stop();
}

void process_collision_events()
{
  simulation::time_event_collision.start();

  #ifdef USE_DEVICE
  #pragma omp target teams distribute parallel for
  #else
  #pragma omp parallel for schedule(runtime)
  #endif
  for (int64_t i = 0; i < simulation::collision_queue.size(); i++) {
    int64_t buffer_idx = simulation::collision_queue[i].idx;
    Particle& p = simulation::device_particles[buffer_idx];
    p.event_collide();
    p.event_revive_from_secondary();
    if (p.alive_)
      dispatch_xs_event(buffer_idx);
  }

  simulation::calculate_fuel_xs_queue.sync_size_device_to_host();
  simulation::calculate_nonfuel_xs_queue.sync_size_device_to_host();
  simulation::collision_queue.resize(0);

  simulation::time_event_collision.stop();
}

void process_death_events(int64_t n_particles)
{
  simulation::time_event_death.start();

  // Local keff tally accumulators
  // (workaround for compilers that don't like reductions w/global variables)
  double absorption = 0.0;
  double collision = 0.0;
  double tracklength = 0.0;
  double leakage = 0.0;

  #ifdef USE_DEVICE
  #pragma omp target teams distribute parallel for reduction(+:absorption, collision, tracklength, leakage)
  #else
  #pragma omp parallel for schedule(runtime) reduction(+:absorption, collision, tracklength, leakage)
  #endif
  for (int64_t i = 0; i < n_particles; i++) {
    Particle& p = simulation::device_particles[i];
    p.accumulate_keff_tallies_local(absorption, collision, tracklength, leakage);
    p.event_death();
  }

  // Write local reduction results to global values
  global_tally_absorption  = absorption;
  global_tally_collision   = collision;
  global_tally_tracklength = tracklength;
  global_tally_leakage     = leakage;

  #ifdef USE_DEVICE

  // Move particle progeny count array back to host
  #pragma omp target update from(simulation::device_progeny_per_particle[:simulation::progeny_per_particle.size()])

  #endif

  simulation::time_event_death.stop();
}

} // namespace openmc
