#include "openmc/event.h"
#include "openmc/material.h"

namespace openmc {

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

void dispatch_xs_event(int64_t i)
{
  Particle* p = &simulation::particles[i];
  int64_t idx;
  if (p->material_ == MATERIAL_VOID) {
    #pragma omp atomic capture
    idx = simulation::calculate_nonfuel_xs_queue_length++;
    simulation::calculate_nonfuel_xs_queue[idx].idx = i;
    simulation::calculate_nonfuel_xs_queue[idx].E = p->E_;
    simulation::calculate_nonfuel_xs_queue[idx].material = p->material_;
    simulation::calculate_nonfuel_xs_queue[idx].type = p->type_;
  }
  else
  {
    if (model::materials[p->material_]->fissionable_) {
      #pragma omp atomic capture
      idx = simulation::calculate_fuel_xs_queue_length++;
      simulation::calculate_fuel_xs_queue[idx].idx = i;
      simulation::calculate_fuel_xs_queue[idx].E = p->E_;
      simulation::calculate_fuel_xs_queue[idx].material = p->material_;
      simulation::calculate_fuel_xs_queue[idx].type = p->type_;
    }
    else
    {
      #pragma omp atomic capture
      idx = simulation::calculate_nonfuel_xs_queue_length++;
      simulation::calculate_nonfuel_xs_queue[idx].idx = i;
      simulation::calculate_nonfuel_xs_queue[idx].E = p->E_;
      simulation::calculate_nonfuel_xs_queue[idx].material = p->material_;
      simulation::calculate_nonfuel_xs_queue[idx].type = p->type_;
    }
  }
}

void process_calculate_xs_events(QueueItem* queue, int64_t n)
{
  // Sort queue by energy
  std::sort(queue, queue+n);

  // Save last_ members, find grid index
  #pragma omp parallel for schedule(runtime)
  for (auto i = 0; i < n; i++) {
    Particle* p = &simulation::particles[queue[i].idx]; 
    p->event_calculate_xs();
  }

  int64_t start = simulation::advance_particle_queue_length;
  int64_t end = start + n;
  int64_t j = 0;
  for (auto i = start; i < end; i++) {
    simulation::advance_particle_queue[i].idx = queue[j].idx;
    simulation::advance_particle_queue[i].E = simulation::particles[queue[j].idx].E_;
    simulation::advance_particle_queue[i].material = simulation::particles[queue[j].idx].material_;
    simulation::advance_particle_queue[i].type = simulation::particles[queue[j].idx].type_;
    j++;
  }
  simulation::advance_particle_queue_length += n;
}

void process_advance_particle_events()
{
  #pragma omp parallel for schedule(runtime)
  for (auto i = 0; i < simulation::advance_particle_queue_length; i++) {
    Particle* p = &simulation::particles[simulation::advance_particle_queue[i].idx];
    p->event_advance();
    if( p->collision_distance_ > p->boundary_.distance ) 
    {
      int64_t idx;
      #pragma omp atomic capture
      idx = simulation::surface_crossing_queue_length++;
      simulation::surface_crossing_queue[idx].idx = simulation::advance_particle_queue[i].idx;
      simulation::surface_crossing_queue[idx].E = p->E_;
      simulation::surface_crossing_queue[idx].material = p->material_;
      simulation::surface_crossing_queue[idx].type = p->type_;
    }
    else
    {
      int64_t idx;
      #pragma omp atomic capture
      idx = simulation::collision_queue_length++;
      simulation::collision_queue[idx].idx = simulation::advance_particle_queue[i].idx;
      simulation::collision_queue[idx].E = p->E_;
      simulation::collision_queue[idx].material = p->material_;
      simulation::collision_queue[idx].type = p->type_;
    }
  }
}

void process_surface_crossing_events()
{
  #pragma omp parallel for schedule(runtime)
  for (auto i = 0; i < simulation::surface_crossing_queue_length; i++) {
    Particle* p = &simulation::particles[simulation::surface_crossing_queue[i].idx];
    p->event_cross_surface();
    p->event_revive_from_secondary();
    if (p->alive_)
      dispatch_xs_event(simulation::surface_crossing_queue[i].idx);
  }

}

void process_collision_events()
{
  #pragma omp parallel for schedule(runtime)
  for (auto i = 0; i < simulation::collision_queue_length; i++) {
    Particle* p = &simulation::particles[simulation::collision_queue[i].idx];
    p->event_collide();
    p->event_revive_from_secondary();
    if (p->alive_)
      dispatch_xs_event(simulation::collision_queue[i].idx);
  }

}

void process_death_events(int64_t n_particles)
{
  #pragma omp parallel for schedule(runtime)
  for (auto i = 0; i < n_particles; i++) {
    Particle& p = simulation::particles[i];
    p.event_death();
  }
}

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

} // namespace openmc
