#ifndef OPENMC_EVENT_H
#define OPENMC_EVENT_H

//! \file event.h
//! \brief Event-based data structures and methods

#include <tuple> // for tie

#include "openmc/particle.h"
#include "openmc/shared_array.h"

namespace openmc {

//==============================================================================
// Structs
//==============================================================================

// In the event-based model, instead of moving or sorting the particles
// themselves based on which event they need, a queue is used to store the
// index (and other useful info) for each event type.
// The EventQueueItem struct holds the relevant information about a particle
// needed for sorting the queue. For very high particle counts, a sorted queue
// has the potential to result in greatly improved cache efficiency. However,
// sorting will introduce some overhead due to the sorting process itself, and
// may not result in any benefits if not enough particles are present for them
// to achieve consistent locality improvements.
struct EventQueueItem {
  int64_t idx;       //!< particle index in event-based particle buffer
  ParticleType type; //!< particle type
  int64_t material;  //!< material that particle is in
  double E;          //!< particle energy

  // Constructors
  EventQueueItem() = default;
  EventQueueItem(const Particle& p, int64_t buffer_idx)
    : idx(buffer_idx), type(p.type()), material(p.material()), E(p.E())
  {}

  // Compare by particle type, then by material type (4.5% fuel/7.0%
  // fuel/cladding/etc), then by energy.
  // TODO: Currently in OpenMC, the material ID corresponds not only to a
  // general type, but also specific isotopic densities. Ideally we would like
  // to be able to just sort by general material type, regardless of densities.
  // A more general material type ID may be added in the future, in which case
  // we can update the material field of this struct to contain the more general
  // id.
  bool operator<(const EventQueueItem& rhs) const
  {
    return std::tie(type, material, E) <
           std::tie(rhs.type, rhs.material, rhs.E);
  }
};

//==============================================================================
// Global variable declarations
//==============================================================================

namespace simulation {

// Event queues. These use the special SharedArray type, rather than a normal
// vector, as they will be shared between threads and may be appended to at the
// same time. To facilitate this, the SharedArray thread_safe_append() method
// is provided which controls the append operations using atomics.
extern SharedArray<EventQueueItem> calculate_fuel_xs_queue;
extern SharedArray<EventQueueItem> calculate_nonfuel_xs_queue;
extern SharedArray<EventQueueItem> advance_particle_queue;
extern SharedArray<EventQueueItem> surface_crossing_queue;
extern SharedArray<EventQueueItem> collision_queue;

// Particle buffer
extern vector<Particle> particles;

} // namespace simulation

//==============================================================================
// Functions
//==============================================================================

//! Allocate space for the event queues and particle buffer
//
//! \param n_particles The number of particles in the particle buffer
void init_event_queues(int64_t n_particles);

//! Free the event queues and particle buffer
void free_event_queues(void);

//! Enqueue a particle based on if it is in fuel or a non-fuel material
//
//! \param buffer_idx The particle's actual index in the particle buffer
void dispatch_xs_event(int64_t buffer_idx);

//! Execute the death event for all particles
//
//! \param n_particles The number of particles in the particle buffer
void process_death_events(int64_t n_particles);

//==============================================================================
// Surface tracking
//==============================================================================

//! Execute the initialization event for all particles
//
//! \param n_particles The number of particles in the particle buffer
//! \param source_offset The offset index in the source bank to use
void process_init_events(int64_t n_particles, int64_t source_offset);

//! Execute the calculate XS event for all particles in this event's buffer
//
//! \param queue A reference to the desired XS lookup queue
void process_calculate_xs_events(SharedArray<EventQueueItem>& queue);

//! Execute the advance particle event for all particles in this event's buffer
void process_advance_particle_events();

//! Execute the surface crossing event for all particles in this event's buffer
void process_surface_crossing_events();

//! Execute the collision event for all particles in this event's buffer
void process_collision_events();

//! Process event queues until all are empty. Each iteration processes the
//! longest queue first to maximize vectorization efficiency.
void process_transport_events();

//! Initialize secondary particles from a shared secondary bank for
//! event-based transport
//
//! \param n_particles The number of particles to initialize
//! \param offset The offset index in the shared secondary bank
//! \param shared_secondary_bank The shared secondary bank to read from
void process_init_secondary_events(int64_t n_particles, int64_t offset,
  const SharedArray<SourceSite>& shared_secondary_bank);

//==============================================================================
// Delta tracking
//==============================================================================

//! Specialization of process_init_events() for delta tracking.
//
//! \param n_particles The number of particles in the particle buffer
//! \param source_offset The offset index in the source bank to use
void process_delta_init_events(int64_t n_particles, int64_t source_offset);

//! Specialization of process_calculate_xs_events() for delta tracking.
//
//! \param queue A reference to the desired XS lookup queue
void process_delta_calculate_xs_events(SharedArray<EventQueueItem>& queue);

//! Execute the delta advance particle event for all particles in this advance
//! buffer
void process_delta_advance_particle_events();

//! Specialization of process_surface_crossing_events() for delta tracking.
void process_delta_surface_crossing_events();

//! Execute the delta tracking collision event for all particles in this event's
//! buffer
void process_delta_collision_events();

//! Specialization of process_transport_events() for delta tracking.
void process_delta_transport_events();

//! Specialization of process_init_secondary_events() for delta tracking.
//
//! \param n_particles The number of particles to initialize
//! \param offset The offset index in the shared secondary bank
//! \param shared_secondary_bank The shared secondary bank to read from
void process_delta_init_secondary_events(int64_t n_particles, int64_t offset,
  const SharedArray<SourceSite>& shared_secondary_bank);

} // namespace openmc

#endif // OPENMC_EVENT_H
