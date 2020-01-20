#ifndef OPENMC_EVENT_H
#define OPENMC_EVENT_H

#include "openmc/particle.h"
#include "openmc/tallies/filter.h"


namespace openmc {

//==============================================================================
// Structs
//==============================================================================

struct QueueItem{
  int64_t idx;         // particle index in event-based buffer
  double E;            // particle energy
  int64_t material;    // material that particle is in
  Particle::Type type; // particle type
  bool operator<(const QueueItem& rhs) const
  {
    // First, compare by particle type
    if (type < rhs.type) {
      return true;
    } else if (type > rhs.type) {
      return false;
    }

    // If we have reached this point, we have the same particle types.
    
    // Next, compare by material:
    // TODO: Currently, material IDs are not usually unique to material
    // types. When unique material type IDs are available, we can alter the
    // material field in this struct to contain the type ID, and then sort
    // by that, as below:
    //if( material < rhs.material)
    //  return true;
    //if( material > rhs.material)
    //  return false;

    // At this point, we have the same particle type, in the same material.
    // Now, compare by energy
    return (E < rhs.E);
  }
};

//==============================================================================
// Global variable declarations
//==============================================================================

namespace simulation {

extern std::unique_ptr<QueueItem[]> calculate_fuel_xs_queue;
extern std::unique_ptr<QueueItem[]> calculate_nonfuel_xs_queue;
extern std::unique_ptr<QueueItem[]> advance_particle_queue;
extern std::unique_ptr<QueueItem[]> surface_crossing_queue;
extern std::unique_ptr<QueueItem[]> collision_queue;
extern std::unique_ptr<Particle[]>  particles;
extern int64_t calculate_fuel_xs_queue_length;
extern int64_t calculate_nonfuel_xs_queue_length;
extern int64_t advance_particle_queue_length;
extern int64_t surface_crossing_queue_length;
extern int64_t collision_queue_length;
extern int64_t max_particles_in_flight;

} // namespace simulation

//==============================================================================
// Functions
//==============================================================================

void init_event_queues(int64_t n_particles);
void free_event_queues(void);
void dispatch_xs_event(int64_t i);
void process_init_events(int64_t n_particles, int64_t source_offset);
void process_calculate_xs_events(QueueItem* queue, int64_t n);
void process_advance_particle_events();
void process_surface_crossing_events();
void process_collision_events();
void process_death_events(int64_t n_particles);

} // namespace openmc

#endif // OPENMC_EVENT_H
