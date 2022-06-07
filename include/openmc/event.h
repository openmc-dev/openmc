#ifndef OPENMC_EVENT_H
#define OPENMC_EVENT_H

//! \file event.h
//! \brief Event-based data structures and methods

#include "openmc/shared_array.h"

namespace openmc {

//==============================================================================
// Structs
//==============================================================================

// Event type macros for queueless mode
#define EVENT_XS_FUEL 0
#define EVENT_XS_NONFUEL 1
#define EVENT_ADVANCE 2
#define EVENT_SURFACE 3
#define EVENT_COLLISION 4
#define EVENT_REVIVAL 5
#define EVENT_DEATH 6

// In the event-based model, instead of moving or sorting the particles
// themselves based on which event they need, a queue is used to store the
// index (and other useful info) for each event type.
// The EventQueueItem struct holds the relevant information about a particle needed
// for sorting the queue. For very high particle counts, a sorted queue has the
// potential to result in greatly improved cache efficiency. However, sorting
// will introduce some overhead due to the sorting process itself, and may not
// result in any benefits if not enough particles are present for them to achieve
// consistent locality improvements. 
struct EventQueueItem{
  int idx;     //!< particle index in event-based particle buffer
  float E;     //!< particle energy
  #ifdef QUEUELESS
  int event;   //!< particle next event type
  #endif

  // Constructors
  EventQueueItem() = default;
  EventQueueItem(double energy, int buffer_idx) :
    idx(buffer_idx), E(static_cast<float>(energy)) {}

  #ifdef COMPILE_CUDA_COMPARATOR
  __host__ __device__
  #endif
  bool operator<(const EventQueueItem& rhs) const
  {
    #ifdef QUEUELESS
    if ( event == rhs.event )
      return E < rhs.E;
    else
      return event < rhs.event;
    #else
    return E < rhs.E;
    #endif
  }
  
  // This is needed by the implementation of parallel quicksort
  bool operator>(const EventQueueItem& rhs) const
  {
    #ifdef QUEUELESS
    if ( event == rhs.event )
      return E > rhs.E;
    else
      return event > rhs.event;
    #else
    return E > rhs.E;
    #endif
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

// Note: we only need to declare the xs queues as global items as the rest
// are only used within the lexical scope of a target construct. 
#pragma omp declare target
#ifdef QUEUELESS
extern SharedArray<EventQueueItem> queue;
#else
extern SharedArray<EventQueueItem> calculate_fuel_xs_queue;
extern SharedArray<EventQueueItem> calculate_nonfuel_xs_queue;
extern SharedArray<EventQueueItem> advance_particle_queue;
extern SharedArray<EventQueueItem> surface_crossing_queue;
extern SharedArray<EventQueueItem> collision_queue;
extern SharedArray<EventQueueItem> revival_queue;
#endif

extern int current_source_offset;
#pragma omp end declare target

extern int sort_counter;

} // namespace simulation

//==============================================================================
// Functions
//==============================================================================

//! Allocate space for the event queues and particle buffer
//
//! \param n_particles The number of particles in the particle buffer
void init_event_queues(int n_particles);

//! Free the event queues and particle buffer
void free_event_queues(void);

//! Enqueue a particle based on if it is in fuel or a non-fuel material
//
//! \param buffer_idx The particle's actual index in the particle buffer
void dispatch_xs_event(int queue_idx, int particle_buffer_idx);

//! Execute the initialization event for all particles
//
//! \param n_particles The number of particles in the particle buffer
void process_init_events(int n_particles);

//! Execute the calculate XS event for all particles in this event's buffer
//
//! \param queue A reference to the desired XS lookup queue
//void process_calculate_xs_events(SharedArray<EventQueueItem>& queue);
void process_calculate_xs_events_fuel(int n_particles);
void process_calculate_xs_events_nonfuel(int n_particles);

//! Execute the advance particle event for all particles in this event's buffer
void process_advance_particle_events(int n_particles);

//! Execute the surface crossing event for all particles in this event's buffer
void process_surface_crossing_events(int n_particles);

//! Execute the collision event for all particles in this event's buffer
void process_collision_events(int n_particles);

//! Execute the death event for all particles
//
//! \param n_particles The number of particles in the particle buffer
void process_death_events(int n_particles);

//! Execute the revival event for all particles in this event's buffer
int process_revival_events(int n_particles, int& n_empty_in_flight_slots);

#ifdef CUDA_THRUST_SORT
//! Sort a queue on-device using CUDA Thrust
//
//! \param begin A pointer to the beginning of the queue
//! \param end A pointer to the end of the queue
void device_sort_event_queue_item(EventQueueItem* begin, EventQueueItem* end);
#endif

#ifdef SYCL_SORT
//! Sort a queue on-device using Intel OneAPI DPL via SYCL interop
//
//! \param begin A pointer to the beginning of the queue
//! \param end A pointer to the end of the queue
void sort_queue_SYCL(EventQueueItem* begin, EventQueueItem* end);
#endif

} // namespace openmc

#endif // OPENMC_EVENT_H
