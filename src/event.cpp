#include "openmc/event.h"
#include "openmc/material.h"

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
std::unique_ptr<Particle[]>  particles;

int64_t calculate_fuel_xs_queue_length {0};
int64_t calculate_nonfuel_xs_queue_length {0};
int64_t advance_particle_queue_length {0};
int64_t surface_crossing_queue_length {0};
int64_t collision_queue_length {0};

int64_t max_particles_in_flight {100000};

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

void enqueue_particle(QueueItem* queue, int64_t& length, Particle* p, int64_t buffer_idx, bool use_atomic)
{
  int64_t idx;
  if (use_atomic) {
    #pragma omp atomic capture
    idx = length++;
  } else {
    idx = length++;
  }

  queue[idx].idx = buffer_idx;
  queue[idx].E = p->E_;
  queue[idx].material = p->material_;
  queue[idx].type = p->type_;
}

void dispatch_xs_event(int64_t i)
{
  Particle* p = &simulation::particles[i];
  int64_t idx;
  if (p->material_ == MATERIAL_VOID ||
      !model::materials[p->material_]->fissionable_) {
    enqueue_particle(simulation::calculate_nonfuel_xs_queue.get(),
        simulation::calculate_nonfuel_xs_queue_length, p, i, true);
  } else {
      enqueue_particle(simulation::calculate_fuel_xs_queue.get(),
          simulation::calculate_fuel_xs_queue_length, p, i, true);
  }
}

void process_calculate_xs_events(QueueItem* queue, int64_t n)
{
  // Sort queue by energy
  std::sort(queue, queue+n);

  // Save last_ members, find grid index
  #pragma omp parallel for schedule(runtime)
  for (int64_t i = 0; i < n; i++) {
    Particle* p = &simulation::particles[queue[i].idx]; 
    p->event_calculate_xs();
    enqueue_particle(simulation::advance_particle_queue.get(),
        simulation::advance_particle_queue_length, p, queue[i].idx, true);
  }
}

void process_advance_particle_events()
{
  #pragma omp parallel for schedule(runtime)
  for (int64_t i = 0; i < simulation::advance_particle_queue_length; i++) {
    int64_t buffer_idx = simulation::advance_particle_queue[i].idx;
    Particle* p = &simulation::particles[buffer_idx];
    p->event_advance();
    if (p->collision_distance_ > p->boundary_.distance) {
      enqueue_particle(simulation::surface_crossing_queue.get(),
          simulation::surface_crossing_queue_length, p, buffer_idx, true);
    } else {
      enqueue_particle(simulation::collision_queue.get(),
          simulation::collision_queue_length, p, buffer_idx, true);
    }
  }
}

void process_surface_crossing_events()
{
  #pragma omp parallel for schedule(runtime)
  for (int64_t i = 0; i < simulation::surface_crossing_queue_length; i++) {
    int64_t buffer_index = simulation::surface_crossing_queue[i].idx;
    Particle* p = &simulation::particles[buffer_index];
    p->event_cross_surface();
    p->event_revive_from_secondary();
    if (p->alive_)
      dispatch_xs_event(buffer_index);
  }

}

void process_collision_events()
{
  #pragma omp parallel for schedule(runtime)
  for (int64_t i = 0; i < simulation::collision_queue_length; i++) {
    int64_t buffer_index = simulation::collision_queue[i].idx;
    Particle* p = &simulation::particles[buffer_index];
    p->event_collide();
    p->event_revive_from_secondary();
    if (p->alive_)
      dispatch_xs_event(buffer_index);
  }

}

void process_death_events(int64_t n_particles)
{
  #pragma omp parallel for schedule(runtime)
  for (int64_t i = 0; i < n_particles; i++) {
    Particle* p = &simulation::particles[i];
    p->event_death();
  }
}

} // namespace openmc
