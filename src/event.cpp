#include "openmc/event.h"
#include "openmc/material.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"
#include "openmc/timer.h"

#ifdef __CUDACC__
#include "openmc/bank.h" // needed to set bank container data in collision kernel
#include "openmc/cuda/advance.h"
#include "openmc/cuda/calculate_xs.h"
#include "openmc/cuda/collide.h"
#include "openmc/cuda/cross_surface.h"
#include "openmc/cuda/death.h"
#include "openmc/cuda/initialize.h"
#include "openmc/cuda/refill.h"
#include "openmc/cuda/util.h" // error handling
#include "openmc/settings.h" // thread_block_size
#include <thrust/sort.h>
#endif

unsigned queue_size;

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
SharedArray<unsigned> dead_particle_indices;

vector<Particle> particles;
pinned_replicated_vector<NuclideMicroXS> micros;

} // namespace simulation

//==============================================================================
// Non-member functions
//==============================================================================

void init_event_queues(unsigned n_particles)
{
  simulation::calculate_fuel_xs_queue.reserve(n_particles);
  simulation::calculate_nonfuel_xs_queue.reserve(n_particles);
  simulation::advance_particle_queue.reserve(n_particles);
  simulation::surface_crossing_queue.reserve(n_particles);
  simulation::collision_queue.reserve(n_particles);

  // If we're not doing SOA particles, allocate an AOS of particles
  // If we are doing SOA particles, those arrays must be allocated
  // after we know how many tallies and nuclides are in the problem.
#ifndef __CUDACC__
  simulation::particles.resize(n_particles);
#endif
  simulation::dead_particle_indices.reserve(n_particles);
  auto tmp = simulation::dead_particle_indices.data();
  cudaMemcpyToSymbol(gpu::dead_particle_indices, &tmp, sizeof(unsigned*));
  queue_size = n_particles;
}

void free_event_queues(void)
{
  simulation::calculate_fuel_xs_queue.clear();
  simulation::calculate_nonfuel_xs_queue.clear();
  simulation::advance_particle_queue.clear();
  simulation::surface_crossing_queue.clear();
  simulation::collision_queue.clear();

  simulation::particles.clear();
  simulation::micros.clear();
}

void dispatch_xs_event(unsigned buffer_idx)
{
  Particle p(buffer_idx);
  if (p.material() == MATERIAL_VOID ||
      !model::materials[p.material()]->fissionable_) {
    simulation::calculate_nonfuel_xs_queue.thread_safe_append({p, buffer_idx});
  } else {
    simulation::calculate_fuel_xs_queue.thread_safe_append({p, buffer_idx});
  }
}

void process_init_events(unsigned n_particles, unsigned source_offset)
{
  simulation::time_event_init.start();

  gpu::managed_calculate_nonfuel_queue_index =
    simulation::calculate_nonfuel_xs_queue.size();
  gpu::managed_calculate_fuel_queue_index =
    simulation::calculate_fuel_xs_queue.size();

  gpu::process_initialize_events_device<BLOCKSIZE>
    <<<n_particles / gpu::thread_block_size + 1, gpu::thread_block_size>>>(
      n_particles, source_offset, simulation::calculate_nonfuel_xs_queue.data(),
      simulation::calculate_fuel_xs_queue.data());
  cudaDeviceSynchronize();
  catchCudaErrors("process_init_events");

  simulation::calculate_nonfuel_xs_queue.updateIndex(
    gpu::managed_calculate_nonfuel_queue_index);
  simulation::calculate_fuel_xs_queue.updateIndex(
    gpu::managed_calculate_fuel_queue_index);

  simulation::time_event_init.stop();
}

void process_calculate_xs_events(SharedArray<EventQueueItem>& queue)
{

  simulation::time_event_sort.start();
  // thrust::sort(queue.begin(), queue.end());
  // cudaDeviceSynchronize();
  simulation::time_event_sort.stop();

  simulation::time_event_calculate_xs.start();
  // TODO: this could possibly be separated into the retrieval of cached
  // XS components and the lookup of new cross sections for nuclides where
  // necessary.
  gpu::process_calculate_xs_events_device<<<
    queue.size() / gpu::thread_block_size + 1, gpu::thread_block_size>>>(
    queue.data(), queue.size());
  cudaDeviceSynchronize();
  catchCudaErrors("process_calculate_xs_events_device");

  auto size_before = simulation::advance_particle_queue.size();
  cudaMemcpy(simulation::advance_particle_queue.end(), queue.begin(),
    queue.size() * sizeof(EventQueueItem), cudaMemcpyDeviceToDevice);
  simulation::advance_particle_queue.updateIndex(size_before + queue.size());
  queue.resize(0);

  simulation::time_event_calculate_xs.stop();
}

void process_advance_particle_events()
{
  simulation::time_event_advance_particle.start();

  // Can't put SharedArrays in managed memory, so these intermediate variables
  // are used to allow pushing back within the kernel. They are both markers
  // for the new size of the queues, and diagnostics in that we'll know if a
  // write out-of-bounds happened after the fact.
  gpu::managed_surface_crossing_queue_index =
    simulation::surface_crossing_queue.size();
  gpu::managed_collision_queue_index = simulation::collision_queue.size();

  gpu::process_advance_events_device<BLOCKSIZE>
    <<<simulation::advance_particle_queue.size() / gpu::thread_block_size + 1,
      gpu::thread_block_size>>>(simulation::advance_particle_queue.data(),
      simulation::advance_particle_queue.size(),
      simulation::surface_crossing_queue.data(),
      simulation::collision_queue.data());
  cudaDeviceSynchronize();
  catchCudaErrors("process_advance_events_device");

  simulation::surface_crossing_queue.updateIndex(
    gpu::managed_surface_crossing_queue_index);
  simulation::collision_queue.updateIndex(gpu::managed_collision_queue_index);

  simulation::advance_particle_queue.resize(0);

  simulation::time_event_advance_particle.stop();
}

void process_surface_crossing_events()
{
  simulation::time_event_surface_crossing.start();

  // Set initial positions of the XS calculation queues for appending
  // while running on GPU
  gpu::managed_calculate_nonfuel_queue_index =
    simulation::calculate_nonfuel_xs_queue.size();
  gpu::managed_calculate_fuel_queue_index =
    simulation::calculate_fuel_xs_queue.size();

  gpu::process_surface_crossing_events_device<BLOCKSIZE>
    <<<simulation::surface_crossing_queue.size() / gpu::thread_block_size + 1,
      gpu::thread_block_size>>>(simulation::surface_crossing_queue.data(),
      simulation::surface_crossing_queue.size(),
      simulation::calculate_nonfuel_xs_queue.data(),
      simulation::calculate_fuel_xs_queue.data());
  cudaDeviceSynchronize();
  catchCudaErrors("process_surface_crossing_events_device");

  simulation::calculate_nonfuel_xs_queue.updateIndex(
    gpu::managed_calculate_nonfuel_queue_index);
  simulation::calculate_fuel_xs_queue.updateIndex(
    gpu::managed_calculate_fuel_queue_index);

  simulation::surface_crossing_queue.resize(0);

  simulation::time_event_surface_crossing.stop();
}

void process_collision_events()
{
  simulation::time_event_collision.start();

  auto fission_bank_start = simulation::fission_bank.data();
  unsigned fission_bank_capacity = simulation::fission_bank.capacity();
  cudaMemcpyToSymbol(
    gpu::fission_bank_start, &fission_bank_start, sizeof(SourceSite*));
  cudaMemcpyToSymbol(
    gpu::fission_bank_capacity, &fission_bank_capacity, sizeof(unsigned));
  gpu::fission_bank_index = simulation::fission_bank.size();

  // Set initial positions of the XS calculation queues for appending
  // while running on GPU
  gpu::managed_calculate_nonfuel_queue_index =
    simulation::calculate_nonfuel_xs_queue.size();
  gpu::managed_calculate_fuel_queue_index =
    simulation::calculate_fuel_xs_queue.size();

  gpu::process_collision_events_device<BLOCKSIZE>
    <<<simulation::collision_queue.size() / gpu::thread_block_size + 1,
      gpu::thread_block_size>>>(simulation::collision_queue.data(),
      simulation::collision_queue.size(),
      simulation::calculate_nonfuel_xs_queue.data(),
      simulation::calculate_fuel_xs_queue.data());
  cudaDeviceSynchronize();
  catchCudaErrors("process_collision_events_device");

  simulation::fission_bank.updateIndex(gpu::fission_bank_index);
  simulation::calculate_nonfuel_xs_queue.updateIndex(
    gpu::managed_calculate_nonfuel_queue_index);
  simulation::calculate_fuel_xs_queue.updateIndex(
    gpu::managed_calculate_fuel_queue_index);

  simulation::collision_queue.resize(0);

  simulation::time_event_collision.stop();
}

unsigned process_refill_events(unsigned remaining_work, unsigned source_offset)
{
  simulation::time_event_refill.start();

  // Firstly, do a compaction on particle indices storing
  // dead particles. This is similar to copy_if, but we want
  // to copy in indices rather than the particles themself.
  simulation::dead_particle_indices.updateIndex(0);
  gpu::dead_particle_indices_indx = 0;
  gpu::scan_for_dead_particles<<<queue_size / gpu::thread_block_size + 1,
    gpu::thread_block_size>>>(queue_size);
  cudaDeviceSynchronize();
  catchCudaErrors("scan_for_dead_particles");
  simulation::dead_particle_indices.updateIndex(
    gpu::dead_particle_indices_indx);

  unsigned num_particles_refilled =
    std::min(simulation::dead_particle_indices.size(), remaining_work);

// Secondly, we loop over dead particle indices, and initialize
// as many fresh particles there as possible.

  gpu::managed_calculate_nonfuel_queue_index =
    simulation::calculate_nonfuel_xs_queue.size();
  gpu::managed_calculate_fuel_queue_index =
    simulation::calculate_fuel_xs_queue.size();
  gpu::refill_dead_particle_slots<BLOCKSIZE>
    <<<num_particles_refilled / gpu::thread_block_size + 1,
      gpu::thread_block_size>>>(num_particles_refilled, source_offset,
      simulation::calculate_nonfuel_xs_queue.data(),
      simulation::calculate_fuel_xs_queue.data());
  cudaDeviceSynchronize();
  catchCudaErrors("refill_dead_particle_slots");
  simulation::calculate_nonfuel_xs_queue.updateIndex(
    gpu::managed_calculate_nonfuel_queue_index);
  simulation::calculate_fuel_xs_queue.updateIndex(
    gpu::managed_calculate_fuel_queue_index);

  simulation::time_event_refill.stop();
  return num_particles_refilled;
}

void process_death_events(unsigned n_particles)
{
  simulation::time_event_death.start();
  // TODO do parallel reduce on particle global tallies here
  gpu::process_death_events_device<<<n_particles / gpu::thread_block_size + 1,
    gpu::thread_block_size>>>(n_particles);
  cudaDeviceSynchronize();
  simulation::time_event_death.stop();
}

} // namespace openmc
