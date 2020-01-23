#ifndef OPENMC_EVENT_H
#define OPENMC_EVENT_H

//! \file event.h
//! \brief Event-based data structures and methods

#include "openmc/particle.h"
#include "openmc/shared_array.h"


namespace openmc {

//==============================================================================
// Structs
//==============================================================================

// In the event-based model, instead of moving or sorting the particles
// themselves based on which event they need, a queue is used to store the
// index (and other useful info) for each event type.
// The QueueItem struct holds the relevant information about a particle needed
// for sorting the queue. For very high particle counts, a sorted queue has the
// potential to result in greatly improved cache efficiency. However, sorting
// will introduce some overhead due to the sorting process itself, and may not
// result in any benefits if not enough particles are present for them to achieve
// consistent locality improvements. 
struct QueueItem{
  int64_t idx;         //!< particle index in event-based particle buffer
  Particle::Type type; //!< particle type
  int64_t material;    //!< material that particle is in
  double E;            //!< particle energy

  // Compare by particle type, then by material type, then by energy
  // TODO: Currently, material IDs are not usually unique to material
  // types. When unique material type IDs are available, we can alter the
  // material field in this struct to contain the material type instead.
  bool operator<(const QueueItem& rhs) const
  {
    return std::tie(type, material, E) < std::tie(rhs.type, rhs.material, rhs.E);
  }
};

//==============================================================================
// Global variable declarations
//==============================================================================

namespace simulation {

// Event queues. These are allocated pointer variables rather than vectors,
// because they are shared between threads and writing to them must be
// coordinated with atomics. This means that normal vector methods (e.g.,
// push_back(), size()) would cause undefined or unintended behavior. Rather,
// adding particles to queues will be done via the enqueue_particle() function.
extern SharedArray<QueueItem> calculate_fuel_xs_queue;
extern SharedArray<QueueItem> calculate_nonfuel_xs_queue;
extern SharedArray<QueueItem> advance_particle_queue;
extern SharedArray<QueueItem> surface_crossing_queue;
extern SharedArray<QueueItem> collision_queue;

// Particle buffer
extern std::vector<Particle>  particles;

} // namespace simulation

//==============================================================================
// Functions
//==============================================================================

//! Allocates space for the event queues and particle buffer
//! \param n_particles The number of particles in the particle buffer
void init_event_queues(int64_t n_particles);

//! Frees the event queues and particle buffer
void free_event_queues(void);

//! Atomically adds a particle to the specified queue
//! \param queue The queue to append the particle to
//! \param p A pointer to the particle
//! \param buffer_idx The particle's actual index in the particle buffer
void enqueue_particle(SharedArray<QueueItem>& queue, const Particle* p, int64_t buffer_idx);

//! Enqueues a particle based on if it is in fuel or a non-fuel material
//! \param buffer_idx The particle's actual index in the particle buffer
void dispatch_xs_event(int64_t buffer_idx);

//! Executes the initialization event for all particles
//! \param n_particles The number of particles in the particle buffer
//! \param source_offset The offset index in the source bank to use
void process_init_events(int64_t n_particles, int64_t source_offset);

//! Executes the calculate XS event for all particles in this event's buffer
//! \param queue A reference to the desired XS lookup queue
void process_calculate_xs_events(SharedArray<QueueItem>& queue);

//! Executes the advance particle event for all particles in this event's buffer
void process_advance_particle_events();

//! Executes the surface crossing event for all particles in this event's buffer
void process_surface_crossing_events();

//! Executes the collision event for all particles in this event's buffer
void process_collision_events();

//! Executes the death event for all particles
//! \param n_particles The number of particles in the particle buffer
void process_death_events(int64_t n_particles);

} // namespace openmc

#endif // OPENMC_EVENT_H
