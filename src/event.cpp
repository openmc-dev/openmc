#include "openmc/bank.h"
#include "openmc/cell.h"
#include "openmc/event.h"
#include "openmc/lattice.h"
#include "openmc/material.h"
#include "openmc/simulation.h"
#include "openmc/surface.h"
#include "openmc/timer.h"
#include "openmc/tallies/tally.h"
#ifndef DEVICE_PRINT
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

std::vector<Particle>  particles;
Particle*  device_particles;

} // namespace simulation

//==============================================================================
// Non-member functions
//==============================================================================

int sort_counter{0};


#pragma omp declare target
int* q_stack {NULL};

void qswap(EventQueueItem& a, EventQueueItem& b)
{
  EventQueueItem tmp = a;
  a = b;
  b = tmp;
}

int partition_by_index(EventQueueItem* arr, int l, int h)
{
  int x = arr[h].idx;
  int i = (l - 1);
  for (int j = l; j <= h - 1; j++) {
    if (arr[j].idx <= x) {
      i++;
      qswap(arr[i], arr[j]);
    }
  }
  qswap(arr[i + 1], arr[h]);
  return (i + 1);
}

int partition_by_energy(EventQueueItem* arr, int l, int h)
{
  double x = arr[h].E;
  int i = (l - 1);
  for (int j = l; j <= h - 1; j++) {
    if (arr[j].E <= x) {
      i++;
      qswap(arr[i], arr[j]);
    }
  }
  qswap(arr[i + 1], arr[h]);
  return (i + 1);
}

void quick_sort(EventQueueItem* arr, int l, int h)
{
  //int q_stack[h-l+1];
  int top = -1;
  q_stack[++top]=l;
  q_stack[++top]=h;
  while(top >= 0)
  {
    h = q_stack[top--];
    l = q_stack[top--];
    int p = partition_by_energy(arr, l, h);
    if(p-1>l)
    {
      q_stack[++top] = l;
      q_stack[++top] = p - 1;
    }
    if(p+1<h)
    {
      q_stack[++top] = p + 1;
      q_stack[++top] = h;
    }
  }
}
#pragma omp end declare target

void quickSort_parallel_internal(EventQueueItem* arr, int left, int right, int cutoff)
{
	int i = left, j = right;
	double tmp;
	double pivot = arr[(left + right) / 2].E;

	{
		while (i <= j) {
			while (arr[i].E < pivot)
				i++;
			while (arr[j].E > pivot)
				j--;
			if (i <= j) {
        qswap(arr[i], arr[j]);
				i++;
				j--;
			}
		}

	}

	if ( ((right-left)<cutoff) ){
		if (left < j){ quickSort_parallel_internal(arr, left, j, cutoff); }
		if (i < right){ quickSort_parallel_internal(arr, i, right, cutoff); }

	}else{
		#pragma omp task
		{ quickSort_parallel_internal(arr, left, j, cutoff); }
		#pragma omp task
		{ quickSort_parallel_internal(arr, i, right, cutoff); }
	}

}

void quickSort_parallel(EventQueueItem* arr, int lenArray){

	// Set minumum problem size to still spawn threads for
	int cutoff = 1000;

	// For this problem size, more than 16 threads on CPU is not helpful
	int	numThreads = 16;

  #pragma omp parallel num_threads(numThreads)
	{
		#pragma omp single nowait
		{
			quickSort_parallel_internal(arr, 0, lenArray-1, cutoff);
		}
	}
}

void sort_queue(SharedArray<EventQueueItem>& queue)
{
  if( simulation::calculate_fuel_xs_queue.size() > settings::sort_frequency )
  {
    {
      #pragma omp target update from(queue.data_[:queue.size()])
      {
        quickSort_parallel(queue.data(), queue.size());
      }
      #pragma omp target update to(queue.data_[:queue.size()])
    }
  }
}

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
  #pragma omp target teams distribute parallel for reduction(+:total_weight)
  #else
  #pragma omp parallel for schedule(runtime) reduction(+:total_weight)
  #endif
  for (int64_t i = 0; i < n_particles; i++) {
    total_weight += simulation::device_particles[i].wgt_;
  }
  simulation::time_event_init.stop();

  // Write total weight to global variable
  simulation::total_weight = total_weight;
  
  simulation::calculate_fuel_xs_queue.sync_size_device_to_host();
  simulation::calculate_nonfuel_xs_queue.sync_size_device_to_host();

}
#ifdef NO_UPDATE
// No Update
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

  int64_t q_size = queue.size();

  #ifdef USE_DEVICE
  #pragma omp target teams distribute parallel for
  #else
  #pragma omp parallel for schedule(runtime)
  #endif
  for (int64_t i = 0; i < q_size; i++) {
    int64_t buffer_idx = queue[i].idx;
    Particle& p = simulation::device_particles[buffer_idx];
    p.event_calculate_xs();
    simulation::advance_particle_queue[offset + i] = queue[i];

    // Updating this on device inside of the kernel saves extra device syncs after launch
    if( i == 0 )
    {
      simulation::advance_particle_queue.size_ = offset + queue.size_;
      queue.size_ = 0;
    }
  }

  // After executing a calculate_xs event, particles will
  // always require an advance event. Therefore, we don't need to use
  // the protected enqueuing function.
  //simulation::advance_particle_queue.resize(offset + queue.size());
  //queue.resize(0);

  simulation::advance_particle_queue.size_ = queue.size_ + offset;
  queue.size_ = 0;
  //#pragma omp target update to(queue.size_, simulation::advance_particle_queue.size_)
  simulation::time_event_calculate_xs.stop();
}
#else

void process_calculate_xs_events_nonfuel()
{
  simulation::time_event_calculate_xs.start();
  
  // The sort here makes less sense, as there are a lot of other various material types, so sorting becomes less powerful
  //sort_queue(simulation::calculate_nonfuel_xs_queue);

  int64_t offset = simulation::advance_particle_queue.size();;

  #ifdef USE_DEVICE
  #pragma omp target teams distribute parallel for
  #else
  #pragma omp parallel for schedule(runtime)
  #endif
  for (int64_t i = 0; i < simulation::calculate_nonfuel_xs_queue.size(); i++) {
    int64_t buffer_idx = simulation::calculate_nonfuel_xs_queue[i].idx;
    Particle& p = simulation::device_particles[buffer_idx];
    p.event_calculate_xs();
    simulation::advance_particle_queue[offset + i] = simulation::calculate_nonfuel_xs_queue[i];
  }

  // After executing a calculate_xs event, particles will
  // always require an advance event. Therefore, we don't need to use
  // the protected enqueuing function.
  simulation::advance_particle_queue.resize(offset + simulation::calculate_nonfuel_xs_queue.size());
  simulation::calculate_nonfuel_xs_queue.resize(0);

  simulation::time_event_calculate_xs.stop();
}

void process_calculate_xs_events_fuel()
{
  simulation::time_event_calculate_xs.start();

  sort_queue(simulation::calculate_fuel_xs_queue);

  int64_t offset = simulation::advance_particle_queue.size();;

  #ifdef USE_DEVICE
  #pragma omp target teams distribute parallel for
  #else
  #pragma omp parallel for schedule(runtime)
  #endif
  for (int64_t i = 0; i < simulation::calculate_fuel_xs_queue.size(); i++) {
    int64_t buffer_idx = simulation::calculate_fuel_xs_queue[i].idx;
    Particle& p = simulation::device_particles[buffer_idx];
    p.event_calculate_xs();
    simulation::advance_particle_queue[offset + i] = simulation::calculate_fuel_xs_queue[i];
  }

  // After executing a calculate_xs event, particles will
  // always require an advance event. Therefore, we don't need to use
  // the protected enqueuing function.
  simulation::advance_particle_queue.resize(offset + simulation::calculate_fuel_xs_queue.size());
  simulation::calculate_fuel_xs_queue.resize(0);

  simulation::time_event_calculate_xs.stop();
}

/*
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
  #pragma omp target teams distribute parallel for
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
*/
#endif

// No update
#ifdef NO_UPDATE
void process_advance_particle_events()
{
  simulation::time_event_advance_particle.start();

  int64_t q_size = simulation::advance_particle_queue.size();

  int hits[2] = {0,0};

  #ifdef USE_DEVICE
  #pragma omp target teams distribute parallel for map(tofrom:hits)
  #else
  #pragma omp parallel for schedule(runtime)
  #endif
  for (int64_t i = 0; i < q_size; i++) {
    int64_t buffer_idx = simulation::advance_particle_queue[i].idx;
    Particle& p = simulation::device_particles[buffer_idx];
    p.event_advance();
    p.event_advance_tally();

    if (p.collision_distance_ > p.boundary_.distance) {
      simulation::surface_crossing_queue.thread_safe_append({p, buffer_idx});
      #pragma omp atomic
      hits[0]++;
    } else {
      simulation::collision_queue.thread_safe_append({p, buffer_idx});
      #pragma omp atomic
      hits[1]++;
    }
    if( i == 0 )
      simulation::advance_particle_queue.size_ = 0;
  }
  //simulation::surface_crossing_queue.sync_size_device_to_host();
  //simulation::collision_queue.sync_size_device_to_host();
  //simulation::advance_particle_queue.resize(0);
  simulation::surface_crossing_queue.size_ += hits[0];
  simulation::collision_queue.size_ += hits[1];
  simulation::advance_particle_queue.size_ = 0;

  simulation::time_event_advance_particle.stop();

}
#else

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

  simulation::time_event_advance_particle.stop();
}
#endif

#ifdef NO_UDPATE
// No update
void process_surface_crossing_events()
{
  simulation::time_event_surface_crossing.start();

  int hits[2] = {0,0};
  int64_t q_size = simulation::surface_crossing_queue.size();

  #ifdef USE_DEVICE
  #pragma omp target teams distribute parallel for map(tofrom:hits)
  #else
  #pragma omp parallel for schedule(runtime)
  #endif
  for (int64_t i = 0; i < q_size; i++) {
    int64_t buffer_idx = simulation::surface_crossing_queue[i].idx;
    Particle& p = simulation::device_particles[buffer_idx];
    p.event_cross_surface();
    p.event_revive_from_secondary();
    if (p.alive_)
    {
      //dispatch_xs_event(buffer_idx);
      if (p.material_ == MATERIAL_VOID || !model::materials[p.material_].fissionable_) {
        simulation::calculate_nonfuel_xs_queue.thread_safe_append({p, buffer_idx});
        #pragma omp atomic
        hits[0]++;
      } else {
        simulation::calculate_fuel_xs_queue.thread_safe_append({p, buffer_idx});
        #pragma omp atomic
        hits[1]++;
      }
    }
    if( i == 0 )
      simulation::surface_crossing_queue.size_ = 0;
  }

  simulation::calculate_nonfuel_xs_queue.size_ += hits[0];
  simulation::calculate_fuel_xs_queue.size_ += hits[1];
  simulation::surface_crossing_queue.size_ = 0;

  //simulation::calculate_fuel_xs_queue.sync_size_device_to_host();
  //simulation::calculate_nonfuel_xs_queue.sync_size_device_to_host();
  //simulation::surface_crossing_queue.resize(0);

  simulation::time_event_surface_crossing.stop();
}
#else

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
#endif

#ifdef NO_UPDATE
// No Update
void process_collision_events()
{
  simulation::time_event_collision.start();

  int hits[2] = {0,0};
  int64_t q_size = simulation::collision_queue.size();

  #ifdef USE_DEVICE
  #pragma omp target teams distribute parallel for map(tofrom:hits)
  #else
  #pragma omp parallel for schedule(runtime)
  #endif
  for (int64_t i = 0; i < q_size; i++) {
    int64_t buffer_idx = simulation::collision_queue[i].idx;
    Particle& p = simulation::device_particles[buffer_idx];
    p.event_collide();
    p.event_revive_from_secondary();
    if (p.alive_)
    {
      //dispatch_xs_event(buffer_idx);
      if (p.material_ == MATERIAL_VOID || !model::materials[p.material_].fissionable_) {
        simulation::calculate_nonfuel_xs_queue.thread_safe_append({p, buffer_idx});
        #pragma omp atomic
        hits[0]++;
      } else {
        simulation::calculate_fuel_xs_queue.thread_safe_append({p, buffer_idx});
        #pragma omp atomic
        hits[1]++;
      }
    }
    if( i == 0 )
      simulation::collision_queue.size_ = 0;
  }

  simulation::calculate_nonfuel_xs_queue.size_ += hits[0];
  simulation::calculate_fuel_xs_queue.size_ += hits[1];
  simulation::collision_queue.size_ = 0;

  //simulation::calculate_fuel_xs_queue.sync_size_device_to_host();
  //simulation::calculate_nonfuel_xs_queue.sync_size_device_to_host();
  //simulation::collision_queue.resize(0);

  simulation::time_event_collision.stop();
}
#else

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
#endif


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
  simulation::time_event_death.stop();

  // Write local reduction results to global values
  global_tally_absorption  = absorption;
  global_tally_collision   = collision;
  global_tally_tracklength = tracklength;
  global_tally_leakage     = leakage;

  #ifdef USE_DEVICE

  // Move particle progeny count array back to host
  #pragma omp target update from(simulation::device_progeny_per_particle[:simulation::progeny_per_particle.size()])

  #endif

}

} // namespace openmc
