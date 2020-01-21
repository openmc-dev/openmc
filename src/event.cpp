#include "openmc/event.h"
#include "openmc/material.h"
#include "openmc/simulation.h"
#include "openmc/timer.h"


namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace simulation {

std::unique_ptr<QueueItem[]> calculate_fuel_xs_queue;
std::unique_ptr<QueueItem[]> calculate_nonfuel_xs_queue;
std::unique_ptr<QueueItem[]> advance_particle_queue;
std::unique_ptr<QueueItem[]> surface_crossing_queue;
std::unique_ptr<QueueItem[]> collision_queue;

int64_t calculate_fuel_xs_queue_length {0};
int64_t calculate_nonfuel_xs_queue_length {0};
int64_t advance_particle_queue_length {0};
int64_t surface_crossing_queue_length {0};
int64_t collision_queue_length {0};

std::unique_ptr<Particle[]>  particles;

} // namespace simulation

//==============================================================================
// Non-member functions
//==============================================================================

void init_event_queues(int64_t n_particles)
{
  simulation::calculate_fuel_xs_queue =    std::make_unique<QueueItem[]>(n_particles);
  simulation::calculate_nonfuel_xs_queue = std::make_unique<QueueItem[]>(n_particles);
  simulation::advance_particle_queue =     std::make_unique<QueueItem[]>(n_particles);
  simulation::surface_crossing_queue =     std::make_unique<QueueItem[]>(n_particles);
  simulation::collision_queue =            std::make_unique<QueueItem[]>(n_particles);
  simulation::particles =                  std::make_unique<Particle[] >(n_particles);
}

void free_event_queues(void)
{
  simulation::calculate_fuel_xs_queue.reset();
  simulation::calculate_nonfuel_xs_queue.reset();
  simulation::advance_particle_queue.reset();
  simulation::surface_crossing_queue.reset();
  simulation::collision_queue.reset();
  simulation::particles.reset();
}

void enqueue_particle(QueueItem* queue, int64_t& length, Particle* p,
    int64_t buffer_idx)
{
  int64_t idx;
  #pragma omp atomic capture
  idx = length++;

  queue[idx].idx = buffer_idx;
  queue[idx].E = p->E_;
  queue[idx].material = p->material_;
  queue[idx].type = p->type_;
}

void dispatch_xs_event(int64_t buffer_idx)
{
  Particle* p = &simulation::particles[buffer_idx];
  if (p->material_ == MATERIAL_VOID ||
      !model::materials[p->material_]->fissionable_) {
    enqueue_particle(simulation::calculate_nonfuel_xs_queue.get(),
        simulation::calculate_nonfuel_xs_queue_length, p, buffer_idx);
  } else {
      enqueue_particle(simulation::calculate_fuel_xs_queue.get(),
          simulation::calculate_fuel_xs_queue_length, p, buffer_idx);
  }
}

void process_init_events(int64_t n_particles, int64_t source_offset)
{
  simulation::time_event_init.start();
  #pragma omp parallel for schedule(runtime)
  for (int64_t i = 0; i < n_particles; i++) {
    initialize_history(&simulation::particles[i], source_offset + i + 1);
    dispatch_xs_event(i);
  }
  simulation::time_event_init.stop();
}

void process_calculate_xs_events(QueueItem* queue, int64_t n_particles)
{
  simulation::time_event_calculate_xs.start();

  // TODO: If using C++17, perform a parallel sort of the queue
  // by particle type, material type, and then energy, in order to
  // improve cache locality and reduce thread divergence on GPU. Prior
  // to C++17, std::sort is a serial only operation, which in this case
  // makes it too slow to be practical for most test problems. 
  //
  // std::sort(std::execution::par_unseq, queue, queue+n);

  #pragma omp parallel for schedule(runtime)
  for (int64_t i = 0; i < n_particles; i++) {
    Particle* p = &simulation::particles[queue[i].idx]; 
    p->event_calculate_xs();
    
    // After executing a calculate_xs event, particles will
    // always require an advance event. Therefore, we don't need to use
    // the protected enqueuing function.
    int64_t offset = simulation::advance_particle_queue_length + i;
    simulation::advance_particle_queue[offset] = queue[i];
  }

  simulation::advance_particle_queue_length += n_particles;

  simulation::time_event_calculate_xs.stop();
}

void process_advance_particle_events()
{
  simulation::time_event_advance_particle.start();

  #pragma omp parallel for schedule(runtime)
  for (int64_t i = 0; i < simulation::advance_particle_queue_length; i++) {
    int64_t buffer_idx = simulation::advance_particle_queue[i].idx;
    Particle* p = &simulation::particles[buffer_idx];
    p->event_advance();
    if (p->collision_distance_ > p->boundary_.distance) {
      enqueue_particle(simulation::surface_crossing_queue.get(),
          simulation::surface_crossing_queue_length, p, buffer_idx);
    } else {
      enqueue_particle(simulation::collision_queue.get(),
          simulation::collision_queue_length, p, buffer_idx);
    }
  }

  simulation::time_event_advance_particle.stop();
}

void process_surface_crossing_events()
{
  simulation::time_event_surface_crossing.start();

  #pragma omp parallel for schedule(runtime)
  for (int64_t i = 0; i < simulation::surface_crossing_queue_length; i++) {
    int64_t buffer_index = simulation::surface_crossing_queue[i].idx;
    Particle* p = &simulation::particles[buffer_index];
    p->event_cross_surface();
    p->event_revive_from_secondary();
    if (p->alive_)
      dispatch_xs_event(buffer_index);
  }
  
  simulation::time_event_surface_crossing.stop();
}

void process_collision_events()
{
  simulation::time_event_collision.start();

  #pragma omp parallel for schedule(runtime)
  for (int64_t i = 0; i < simulation::collision_queue_length; i++) {
    int64_t buffer_index = simulation::collision_queue[i].idx;
    Particle* p = &simulation::particles[buffer_index];
    p->event_collide();
    p->event_revive_from_secondary();
    if (p->alive_)
      dispatch_xs_event(buffer_index);
  }

  simulation::time_event_collision.stop();
}

void process_death_events(int64_t n_particles)
{
  simulation::time_event_death.start();
  #pragma omp parallel for schedule(runtime)
  for (int64_t i = 0; i < n_particles; i++) {
    Particle* p = &simulation::particles[i];
    p->event_death();
  }
  simulation::time_event_death.stop();
}

} // namespace openmc
