#include "openmc/simulation.h"

#include "openmc/bank.h"
#include "openmc/capi.h"
#include "openmc/cell.h"
#include "openmc/container_util.h"
#include "openmc/eigenvalue.h"
#include "openmc/error.h"
#include "openmc/geometry.h"
#include "openmc/material.h"
#include "openmc/message_passing.h"
#include "openmc/mgxs_interface.h"
#include "openmc/nuclide.h"
#include "openmc/output.h"
#include "openmc/particle.h"
#include "openmc/photon.h"
#include "openmc/physics.h"
#include "openmc/physics_mg.h"
#include "openmc/random_lcg.h"
#include "openmc/settings.h"
#include "openmc/source.h"
#include "openmc/state_point.h"
#include "openmc/thermal.h"
#include "openmc/timer.h"
#include "openmc/tallies/derivative.h"
#include "openmc/tallies/filter.h"
#include "openmc/tallies/tally.h"
#include "openmc/tallies/tally_scoring.h"
#include "openmc/tallies/trigger.h"
#include "openmc/track_output.h"

#ifdef _OPENMP
#include <omp.h>
#endif
#include "xtensor/xview.hpp"

#ifdef OPENMC_MPI
#include <mpi.h>
#endif

#include <algorithm>
#include <string>
#include <chrono> 

#define DYNAMIC_SIZE 8

namespace openmc {

	/*
extern std::vector<Particle*> calculate_fuel_xs_queue;
extern std::vector<Particle*> calculate_nonfuel_xs_queue;
extern std::vector<Particle*> advance_particle_queue;
extern std::vector<Particle*> surface_crossing_queue;
extern std::vector<Particle*> collision_queue;
#pragma omp threadprivate(calculate_fuel_xs_queue, calculate_nonfuel_xs_queue, advance_particle_queue, surface_crossing_queue, collision_queue)

std::vector<Particle*> calculate_fuel_xs_queue;
std::vector<Particle*> calculate_nonfuel_xs_queue;
std::vector<Particle*> advance_particle_queue;
std::vector<Particle*> surface_crossing_queue;
std::vector<Particle*> collision_queue;
*/

struct QueueItem{
	int idx;      // particle index in event-based buffer
	double E;     // particle energy
	int material; // material that particle is in
  Particle::Type type;
  bool operator<(const QueueItem & rhs) const
  {
    // First, compare by type
    if( type < rhs.type )
      return true;
    if( type > rhs.type )
      return false;

    // At this point, we have the same particle types.
    // Now, compare by material
    
    // TODO: Temporarily disabled as SMR problem has different material IDs for every pin
    // Need to sort by material type instead...
    /*
    if( material < rhs.material)
      return true;
    if( material > rhs.material)
      return false;
      */

    // At this point, we have the same particle type, in the same material.
    // Now, compare by energy
    return (E < rhs.E);
  }
};
bool by_energy   (QueueItem a, QueueItem b) { return (a.E        < b.E); }
bool by_material (QueueItem a, QueueItem b) { return (a.material < b.material); }

QueueItem * calculate_fuel_xs_queue;
QueueItem * calculate_nonfuel_xs_queue;
QueueItem * advance_particle_queue;
QueueItem * surface_crossing_queue;
QueueItem * collision_queue;
Particle * particles;

int calculate_fuel_xs_queue_length    = 0;
int calculate_nonfuel_xs_queue_length = 0;
int advance_particle_queue_length     = 0;
int surface_crossing_queue_length     = 0;
int collision_queue_length            = 0;

const int MAX_PARTICLES_IN_FLIGHT = 1000000;

void init_event_queues(int n_particles)
{
	calculate_fuel_xs_queue =    new QueueItem[n_particles];
	calculate_nonfuel_xs_queue = new QueueItem[n_particles];
	advance_particle_queue =     new QueueItem[n_particles];
	surface_crossing_queue =     new QueueItem[n_particles];
	collision_queue =            new QueueItem[n_particles];
	particles =                  new Particle[n_particles];
}

void free_event_queues(void)
{
	delete[] calculate_fuel_xs_queue;
	delete[] calculate_nonfuel_xs_queue;
	delete[] advance_particle_queue;
	delete[] surface_crossing_queue;
	delete[] collision_queue;
	delete[] particles;
}

constexpr size_t MAX_PARTICLES_PER_THREAD {100};

Particle::Bank * shared_fission_bank;
int shared_fission_bank_length = 0;
int shared_fission_bank_max;

void init_shared_fission_bank(int max)
{
	shared_fission_bank_max = max;
	shared_fission_bank = new Particle::Bank[max];
}

void free_shared_fission_bank(void)
{
	delete[] shared_fission_bank;
	shared_fission_bank_length = 0;
}

// TODO: What is going on here?
void revive_particle_from_secondary(Particle* p)
{
  //p->from_source(&simulation::secondary_bank.back());
  p->from_source(&p->secondary_bank_.back());
  //simulation::secondary_bank.pop_back();
  p->secondary_bank_.pop_back();
  // n_event = 0;

  // Enter new particle in particle track file
  if (p->write_track_) add_particle_track(*p);
}

void dispatch_xs_event(int i)
{
  Particle * p = particles + i;
  int idx;
  if (p->material_ == MATERIAL_VOID) {
    #pragma omp atomic capture
    idx = calculate_nonfuel_xs_queue_length++;
    //std::cout << "Dispatching particle to non Fuel XS queue idx = " << idx << std::endl;
    calculate_nonfuel_xs_queue[idx].idx = i;
    calculate_nonfuel_xs_queue[idx].E = p->E_;
    calculate_nonfuel_xs_queue[idx].material = p->material_;
    calculate_nonfuel_xs_queue[idx].type = p->type_;
  }
  else
  {
    if (model::materials[p->material_]->fissionable_) {
      #pragma omp atomic capture
      idx = calculate_fuel_xs_queue_length++;
      //std::cout << "Dispatching particle to Fuel XS queue idx = " << idx << std::endl;
      calculate_fuel_xs_queue[idx].idx = i;
      calculate_fuel_xs_queue[idx].E = p->E_;
      calculate_fuel_xs_queue[idx].material = p->material_;
      calculate_fuel_xs_queue[idx].type = p->type_;
    }
    else
    {
      #pragma omp atomic capture
      idx = calculate_nonfuel_xs_queue_length++;
      //std::cout << "Dispatching particle to non Fuel XS queue idx = " << idx << std::endl;
      calculate_nonfuel_xs_queue[idx].idx = i;
      calculate_nonfuel_xs_queue[idx].E = p->E_;
      calculate_nonfuel_xs_queue[idx].material = p->material_;
      calculate_nonfuel_xs_queue[idx].type = p->type_;
    }
  }
}

void process_calculate_xs_events(QueueItem * queue, int n)
{
  // Sort queue by energy
  std::sort(queue, queue+n);
  //std::sort(queue, queue+n, by_energy);

  // Then, stable sort by material (so as to preserve energy ordering)
  //std::stable_sort(queue, queue+n, by_material);

  // Save last_ members, find grid index
  int lost_particles = 0;
  #pragma omp parallel for reduction(+:lost_particles) schedule(dynamic,DYNAMIC_SIZE)
  for (int i = 0; i < n; i++) {
	  Particle *p = particles + queue[i].idx; 
	  //std::cout << "particle offset = " << queue[i] << std::endl;
    // Set the random number stream
	// TODO: Move RNG seeds to particle storage
  /*
    if (p->type_ == Particle::Type::neutron) {
		p->stream_ = STREAM_TRACKING;
    } else {
		p->stream_ = STREAM_PHOTON;
    }

    // Store pre-collision particle properties
    p->wgt_last_ = p->wgt_;
    p->E_last_ = p->E_;
    p->u_last_ = p->u();
    p->r_last_ = p->r();
    
    // Reset event variables
    p->event_ = EVENT_KILL;
    p->event_nuclide_ = NUCLIDE_NONE;
    p->event_mt_ = REACTION_NONE;

    // If the cell hasn't been determined based on the particle's location,
    // initiate a search for the current cell. This generally happens at the
    // beginning of the history and again for any secondary particles
    if (p->coord_[p->n_coord_ - 1].cell == C_NONE) {
      if (!find_cell(p, false)) {
        p->mark_as_lost("Could not find the cell containing particle "
          + std::to_string(p->id_));
        //return;
		lost_particles += 1;
		continue;
      }

      // set birth cell attribute
      if (p->cell_born_ == C_NONE) p->cell_born_ = p->coord_[p->n_coord_ - 1].cell;
    }

    // Write particle track.
    if (p->write_track_) write_particle_track(*p);

    if (settings::check_overlaps) check_cell_overlap(p);
    */
    p->event_calculate_xs_I();
  }

  if( lost_particles > 0 )
	  exit(1);

  #pragma omp parallel for schedule(dynamic, DYNAMIC_SIZE)
  for( int i = 0; i < n; i++ )
  {
	  Particle * p = particles + queue[i].idx;
    /*
	  // Calculate microscopic and macroscopic cross sections
	  if (p->material_ != MATERIAL_VOID) {
		  if (settings::run_CE) {
			  if (p->material_ != p->material_last_ || p->sqrtkT_ != p->sqrtkT_last_) {
				  // If the material is the same as the last material and the
				  // temperature hasn't changed, we don't need to lookup cross
				  // sections again.
				  model::materials[p->material_]->calculate_xs(*p);
			  }
		  } // else MG not supported
	  } else {
		  p->macro_xs_.total      = 0.0;
		  p->macro_xs_.absorption = 0.0;
		  p->macro_xs_.fission    = 0.0;
		  p->macro_xs_.nu_fission = 0.0;
	  }
    */
    p->event_calculate_xs_II();
  }

  int start = advance_particle_queue_length;
  int end = start + n;
  int j = 0;
  for( int i = start; i < end; i++ )
  {
	  advance_particle_queue[i].idx = queue[j].idx;
	  advance_particle_queue[i].E = particles[queue[j].idx].E_;
	  advance_particle_queue[i].material = particles[queue[j].idx].material_;
	  advance_particle_queue[i].type = particles[queue[j].idx].type_;
	  j++;
  }
  advance_particle_queue_length += n;
}

void process_advance_particle_events()
{
  //for (auto& p : advance_particle_queue) {
  #pragma omp parallel for schedule(dynamic, DYNAMIC_SIZE)
  for (int i = 0; i < advance_particle_queue_length; i++) {
	  Particle * p = particles + advance_particle_queue[i].idx;
    //p->trace_ == (p->id_ == 0);

    /*
    // Sample a distance to collision
    double d_collision;
    if (p->type_ == Particle::Type::electron ||
        p->type_ == Particle::Type::positron) {
      d_collision = 0.0;
    } else if (p->macro_xs_.total == 0.0) {
      d_collision = INFINITY;
    } else {
      d_collision = -std::log(prn(p->current_seed())) / p->macro_xs_.total;
    }

    // -------------- break here? -------------------

    // Find the distance to the nearest boundary
    p->boundary_ = distance_to_boundary(p);

    // Select smaller of the two distances
    double distance;
	int idx;
    if (p->boundary_.distance < d_collision) {
		#pragma omp atomic capture
		idx = surface_crossing_queue_length++;
      surface_crossing_queue[idx].idx = advance_particle_queue[i].idx;
      surface_crossing_queue[idx].E = p->E_;
      surface_crossing_queue[idx].material = p->material_;
      surface_crossing_queue[idx].type = p->type_;
      distance = p->boundary_.distance;
    } else {
		#pragma omp atomic capture
		idx = collision_queue_length++;
		collision_queue[idx].idx = advance_particle_queue[i].idx;
		collision_queue[idx].E = p->E_;
		collision_queue[idx].material = p->material_;
		collision_queue[idx].type = p->type_;
      distance = d_collision;
    }

    // -------------- break here? -------------------

    // Advance particle
    for (int j = 0; j < p->n_coord_; ++j) {
      p->coord_[j].r += distance * p->coord_[j].u;
    }

    // -------------- break here? -------------------

    // Score track-length tallies
    if (!model::active_tracklength_tallies.empty()) {
      score_tracklength_tally(p, distance);
    }

    // Score track-length estimate of k-eff
    if (settings::run_mode == RUN_MODE_EIGENVALUE &&
        p->type_ == Particle::Type::neutron) {
      p->tally_tracklength_ += p->wgt_ * distance * p->macro_xs_.nu_fission;
    }

    // Score flux derivative accumulators for differential tallies.
    if (!model::active_tallies.empty()) {
      score_track_derivative(p, distance);
    }
  */
    p->event_advance();
    if( p->collision_distance_ > p->boundary_.distance ) 
    {
      int idx;
      #pragma omp atomic capture
      idx = surface_crossing_queue_length++;
      surface_crossing_queue[idx].idx = advance_particle_queue[i].idx;
      surface_crossing_queue[idx].E = p->E_;
      surface_crossing_queue[idx].material = p->material_;
      surface_crossing_queue[idx].type = p->type_;
    }
    else
    {
      int idx;
      #pragma omp atomic capture
      idx = collision_queue_length++;
      collision_queue[idx].idx = advance_particle_queue[i].idx;
      collision_queue[idx].E = p->E_;
      collision_queue[idx].material = p->material_;
      collision_queue[idx].type = p->type_;
    }
  }
  advance_particle_queue_length = 0;
}

void process_surface_crossing_events()
{
  //for (auto& p : surface_crossing_queue) {
  #pragma omp parallel for schedule(dynamic, DYNAMIC_SIZE)
  for (int i = 0; i < surface_crossing_queue_length; i++) {
	  Particle * p = particles + surface_crossing_queue[i].idx;
    /*
    // Set surface that particle is on and adjust coordinate levels
    p->surface_ = p->boundary_.surface_index;
    p->n_coord_ = p->boundary_.coord_level;

    // Saving previous cell data
    for (int j = 0; j < p->n_coord_; ++j) {
      p->cell_last_[j] = p->coord_[j].cell;
    }
    p->n_coord_last_ = p->n_coord_;

    if (p->boundary_.lattice_translation[0] != 0 ||
        p->boundary_.lattice_translation[1] != 0 ||
        p->boundary_.lattice_translation[2] != 0) {
      // Particle crosses lattice boundary
      cross_lattice(p, p->boundary_);
      p->event_ = EVENT_LATTICE;
    } else {
      // Particle crosses surface
      p->cross_surface();
      p->event_ = EVENT_SURFACE;
    }
    // Score cell to cell partial currents
    if (!model::active_surface_tallies.empty()) {
      score_surface_tally(p, model::active_surface_tallies);
    }

    //if (!p->alive_ && !simulation::secondary_bank.empty()) {
    if (!p->alive_ && !p->secondary_bank_.empty()) {
      revive_particle_from_secondary(p);
    }
    */
    p->event_cross_surface();

    if (p->alive_)
    {
      dispatch_xs_event(surface_crossing_queue[i].idx);
    }
  }

  surface_crossing_queue_length = 0;
}

void process_collision_events()
{
  //for (auto& p : collision_queue) {
  #pragma omp parallel for schedule(dynamic,DYNAMIC_SIZE)
  for (int i = 0; i < collision_queue_length; i++) {
	  Particle * p = particles + collision_queue[i].idx;
    /*
	  //std::cout << "Beginning collision of particle id " << collision_queue[i] << " with energy E = " << p->E_ << std::endl;
    // Score collision estimate of keff
    if (settings::run_mode == RUN_MODE_EIGENVALUE &&
        p->type_ == Particle::Type::neutron) {
      p->tally_collision_ += p->wgt_ * p->macro_xs_.nu_fission
        / p->macro_xs_.total;
    }

    // Score surface current tallies -- this has to be done before the collision
    // since the direction of the particle will change and we need to use the
    // pre-collision direction to figure out what mesh surfaces were crossed

    if (!model::active_meshsurf_tallies.empty())
      score_surface_tally(p, model::active_meshsurf_tallies);
	    
	//std::cout << "After surface tally of particle id " << collision_queue[i] << " with energy E = " << p->E_ << std::endl;

    // Clear surface component
    p->surface_ = 0;

    if (settings::run_CE) {
      collision(p);
    } else {
      collision_mg(p);
    }

	//std::cout << "After collision() of particle id " << collision_queue[i] << " with energy E = " << p->E_ << std::endl;

    // Score collision estimator tallies -- this is done after a collision
    // has occurred rather than before because we need information on the
    // outgoing energy for any tallies with an outgoing energy filter
    if (!model::active_collision_tallies.empty()) score_collision_tally(p);
    if (!model::active_analog_tallies.empty()) {
      if (settings::run_CE) {
        score_analog_tally_ce(p);
      } else {
        score_analog_tally_mg(p);
      }
    }
	//std::cout << "After analog tally of particle id " << collision_queue[i] << " with energy E = " << p->E_ << std::endl;

    // Reset banked weight during collision
    p->n_bank_ = 0;
    p->n_bank_second_ = 0;
    p->wgt_bank_ = 0.0;
    for (int& v : p->n_delayed_bank_) v = 0;

    // Reset fission logical
    p->fission_ = false;

    // Save coordinates for tallying purposes
    p->r_last_current_ = p->r();

    // Set last material to none since cross sections will need to be
    // re-evaluated
    p->material_last_ = C_NONE;

    // Set all directions to base level -- right now, after a collision, only
    // the base level directions are changed
    for (int j = 0; j < p->n_coord_ - 1; ++j) {
      if (p->coord_[j + 1].rotated) {
        // If next level is rotated, apply rotation matrix
        const auto& m {model::cells[p->coord_[j].cell]->rotation_};
        const auto& u {p->coord_[j].u};
        p->coord_[j + 1].u = u.rotate(m);
      } else {
        // Otherwise, copy this level's direction
        p->coord_[j+1].u = p->coord_[j].u;
      }
    }

    // Score flux derivative accumulators for differential tallies.
    if (!model::active_tallies.empty()) score_collision_derivative(p);
    */
    p->event_collide();

    //if (!p->alive_ && !simulation::secondary_bank.empty()) {
    /*
    if (!p->alive_ && !p->secondary_bank_.empty()) {
      revive_particle_from_secondary(p);
    }
    */

    p->event_revive_from_secondary();

    if (p->alive_)
    {
      dispatch_xs_event(collision_queue[i].idx);
      //std::cout << "Ended collision of particle id " << collision_queue[i] << " with energy E = " << p->E_ << std::endl;
      //assert(std::isfinite(p->E_) );
    }
  }

  collision_queue_length = 0;
}

/*
void check_energies(void)
{
	int * Q;
	int n;


	Q = calculate_fuel_xs_queue;
	n = calculate_fuel_xs_queue_length;
	for( int i = 0; i < n; i++ )
	{
		if( !std::isfinite(particles[Q[i]].E_ ) )
		{
			std::cout << "NAN energy particle found at index xs FUEL " << Q[i] << std::endl;
			assert(0);
		}
	}
	Q = calculate_nonfuel_xs_queue;
	n = calculate_nonfuel_xs_queue_length ;
	for( int i = 0; i < n; i++ )
	{
		if( !std::isfinite(particles[Q[i]].E_ ) )
		{
			std::cout << "NAN energy particle found at index xs Non fuel " << Q[i] << std::endl;
			assert(0);
		}
	}
	Q = advance_particle_queue;
	n = advance_particle_queue_length     ;
	for( int i = 0; i < n; i++ )
	{
		if( !std::isfinite(particles[Q[i]].E_ ) )
		{
			std::cout << "NAN energy particle found at index advance particle " << Q[i] << std::endl;
			assert(0);
		}
	}
	Q = surface_crossing_queue;
	n = surface_crossing_queue_length     ;
	for( int i = 0; i < n; i++ )
	{
		if( !std::isfinite(particles[Q[i]].E_ ) )
		{
			std::cout << "NAN energy particle found at index surface crossing " << Q[i] << std::endl;
			assert(0);
		}
	}
	Q = collision_queue;
	n = collision_queue_length            ;
	for( int i = 0; i < n; i++ )
	{
		if( !std::isfinite(particles[Q[i]].E_ ) )
		{
			std::cout << "NAN energy particle found at index collision " << Q[i] << std::endl;
			assert(0);
		}
	}
}
*/

double get_time()
{
	#ifdef _OPENMP
	return omp_get_wtime();
	#endif

	#ifdef OPENMC_MPI
	return MPI_Wtime();
	#endif

	unsigned long us_since_epoch = std::chrono::high_resolution_clock::now().time_since_epoch() / std::chrono::microseconds(1);
	return (double) us_since_epoch / 1.0e6;
}

void transport()
{
	double stop, start;
  double time_init = 0;
	double time_fuel_xs = 0;
	double time_nonfuel_xs = 0;
	double time_advance = 0;
	double time_collision = 0;
	double time_surf = 0;


  start = get_time();

	int remaining_work = simulation::work_per_rank;
	int source_offset = 0;
		
	int max_n_particles = MAX_PARTICLES_IN_FLIGHT;
	if( max_n_particles > remaining_work)
		max_n_particles = remaining_work;
	init_event_queues(max_n_particles);

  stop = get_time();
  time_init += stop - start;

	// Subiterations to complete sets of particles
	while (remaining_work > 0) {
    start = get_time();


		// Figure out work for this subiteration
		int n_particles = MAX_PARTICLES_IN_FLIGHT;
		if( n_particles > remaining_work)
			n_particles = remaining_work;

		//std::cout << "Initializing particle histories..." << std::endl;
		// Initialize all histories
		// TODO: Parallelize

		/*
		for( int p = 0; p < mpi::n_procs; p++)
		{
			MPI_Barrier(mpi::intracomm);
			if( p == mpi::rank )
			{
			*/
        #pragma omp parallel for schedule(dynamic, DYNAMIC_SIZE)
				for (int i = 0; i < n_particles; i++) {
					initialize_history(particles + i, source_offset + i + 1);
				}
				/*
			}
		}
		*/

		//std::cout << "Enqueing particles for XS Lookups..." << std::endl;
		// Add all particles to advance particle queue
		// TODO: Parallelize
     #pragma omp parallel for schedule(dynamic, DYNAMIC_SIZE)
		for (int i = 0; i < n_particles; i++) {
			dispatch_xs_event(i);
		}
  
    stop = get_time();
    time_init += stop - start;

		int event_kernel_executions = 0;
		while (true) {
			event_kernel_executions++;
      /*
			std::cout << "Fuel XS Lookups = " << calculate_fuel_xs_queue_length << std::endl;
			std::cout << "Non Fuel XS Lookups = " << calculate_nonfuel_xs_queue_length << std::endl;
			std::cout << "Advance Particles = " << advance_particle_queue_length << std::endl;
			std::cout << "Surface Crossings = " << surface_crossing_queue_length << std::endl;
			std::cout << "Collisions = " << collision_queue_length << std::endl;
      */
			/*
			Particle * p = particles +  1;
		   std::cout << "E = " << p->E_ << " and Position {" <<
			   p->r().x << ", " <<
			   p->r().y << ", " <<
			   p->r().z << "}" << std::endl;
			   */
			int max = std::max({calculate_fuel_xs_queue_length, calculate_nonfuel_xs_queue_length, advance_particle_queue_length, surface_crossing_queue_length, collision_queue_length});
			//check_energies();
			if (max == 0) {
				break;
			} else if (max == calculate_fuel_xs_queue_length) {
				//std::cout << "pre fuel XS check..." << std::endl;
				//check_energies(calculate_fuel_xs_queue, calculate_fuel_xs_queue_length);
				//std::cout << "Performing Fuel XS Lookups..." << std::endl;
				start = get_time();
				process_calculate_xs_events(calculate_fuel_xs_queue, calculate_fuel_xs_queue_length);
				stop = get_time();
				time_fuel_xs += (stop-start);
				calculate_fuel_xs_queue_length = 0;
			} else if (max == calculate_nonfuel_xs_queue_length) {
				//std::cout << "pre non fuel XS check..." << std::endl;
				//check_energies(calculate_nonfuel_xs_queue, calculate_nonfuel_xs_queue_length);
			//	std::cout << "Performing Non Fuel XS Lookups..." << std::endl;
				start = get_time();
				process_calculate_xs_events(calculate_nonfuel_xs_queue, calculate_nonfuel_xs_queue_length);
				stop = get_time();
				time_nonfuel_xs += (stop-start);
				calculate_nonfuel_xs_queue_length = 0;
			} else if (max == advance_particle_queue_length) {
				//std::cout << "pre advancing check..." << std::endl;
				//check_energies(advance_particle_queue, advance_particle_queue_length);
				//std::cout << "Advancing Particles..." << std::endl;
				start = get_time();
				process_advance_particle_events();
				stop = get_time();
				time_advance += (stop-start);
			} else if (max == surface_crossing_queue_length) {
				//std::cout << "pre surface crossing check..." << std::endl;
				//check_energies(surface_crossing_queue, surface_crossing_queue_length);
				//std::cout << "Surface Crossings..." << std::endl;
				start = get_time();
				process_surface_crossing_events();
				stop = get_time();
				time_surf += (stop-start);
			} else if (max == collision_queue_length) {
				//std::cout << "pre Colliding check..." << std::endl;
				//check_energies(collision_queue, collision_queue_length);
				//std::cout << "Colliding..." << std::endl;
				start = get_time();
				process_collision_events();
				stop = get_time();
				time_collision += (stop-start);
			}
		}

    // Finish particle track output and contribute to global tally variables
    for (int i = 0; i < n_particles; i++) {
      Particle& p = particles[i];
      if (p.write_track_) {
        write_particle_track(p);
        finalize_particle_track(p);
      }
      global_tally_absorption +=  p.tally_absorption_;
      global_tally_collision +=   p.tally_collision_;
      global_tally_tracklength += p.tally_tracklength_;
      global_tally_leakage +=     p.tally_leakage_;
    }

		remaining_work -= n_particles;
		source_offset += n_particles;

		// Should all be zero
		/*
		calculate_fuel_xs_queue_length    = 0;
		calculate_nonfuel_xs_queue_length = 0;
		advance_particle_queue_length     = 0;
		surface_crossing_queue_length     = 0;
		collision_queue_length            = 0;
		*/
		
		std::cout << "Event kernels retired: " << event_kernel_executions << std::endl;
	}
	if( mpi::rank == 0 )
	{
		std::cout << "Particle Init Time: " << time_init << std::endl;
		std::cout << "Fuel XS Time:       " << time_fuel_xs << std::endl;
		std::cout << "Non Fuel XS Time:   " << time_nonfuel_xs << std::endl;
		std::cout << "Advance Time:       " << time_advance << std::endl;
		std::cout << "Surface Time:       " << time_surf << std::endl;
		std::cout << "Collision Time:     " << time_collision<< std::endl;
	}
	//shared_fission_bank_length = 0;
	free_event_queues();
}

} // namespace openmc

//==============================================================================
// C API functions
//==============================================================================

// OPENMC_RUN encompasses all the main logic where iterations are performed
// over the batches, generations, and histories in a fixed source or k-eigenvalue
// calculation.

int openmc_run()
{
  openmc_simulation_init();

  int err = 0;
  int status = 0;
  while (status == 0 && err == 0) {
    err = openmc_next_batch(&status);
  }

  openmc_simulation_finalize();
  return err;
}

int openmc_simulation_init()
{
  using namespace openmc;

  // Skip if simulation has already been initialized
  if (simulation::initialized) return 0;

  // Determine how much work each process should do
  calculate_work();

  // Allocate source bank, and for eigenvalue simulations also allocate the
  // fission bank
  allocate_banks();
  init_shared_fission_bank(simulation::work_per_rank * 3);

  // Allocate tally results arrays if they're not allocated yet

  for (auto& t : model::tallies) {
    t->init_results();
  }

  // Set up material nuclide index mapping
  for (auto& mat : model::materials) {
    mat->init_nuclide_index();
  }

  // Reset global variables -- this is done before loading state point (as that
  // will potentially populate k_generation and entropy)
  simulation::current_batch = 0;
  simulation::k_generation.clear();
  simulation::entropy.clear();
  simulation::need_depletion_rx = false;
  openmc_reset();

  // If this is a restart run, load the state point data and binary source
  // file
  if (settings::restart_run) {
    load_state_point();
    write_message("Resuming simulation...", 6);
  } else {
    initialize_source();
  }

  // Display header
  if (mpi::master) {
    if (settings::run_mode == RUN_MODE_FIXEDSOURCE) {
      header("FIXED SOURCE TRANSPORT SIMULATION", 3);
    } else if (settings::run_mode == RUN_MODE_EIGENVALUE) {
      header("K EIGENVALUE SIMULATION", 3);
      if (settings::verbosity >= 7) print_columns();
    }
  }

  // Set flag indicating initialization is done
  simulation::initialized = true;
  return 0;
}

int openmc_simulation_finalize()
{
  using namespace openmc;

  // Skip if simulation was never run
  if (!simulation::initialized) return 0;

  // Stop active batch timer and start finalization timer
  simulation::time_active.stop();
  simulation::time_finalize.start();

  // Clear material nuclide mapping
  for (auto& mat : model::materials) {
    mat->mat_nuclide_index_.clear();
  }

  // Increment total number of generations
  simulation::total_gen += simulation::current_batch*settings::gen_per_batch;

#ifdef OPENMC_MPI
  broadcast_results();
#endif

  // Write tally results to tallies.out
  if (settings::output_tallies && mpi::master) write_tallies();

  // Deactivate all tallies
  for (auto& t : model::tallies) {
    t->active_ = false;
  }

  // Stop timers and show timing statistics
  simulation::time_finalize.stop();
  simulation::time_total.stop();
  if (mpi::master) {
    if (settings::verbosity >= 6) print_runtime();
    if (settings::verbosity >= 4) print_results();
  }
  if (settings::check_overlaps) print_overlap_check();

  free_shared_fission_bank();

  // Reset flags
  simulation::need_depletion_rx = false;
  simulation::initialized = false;
  return 0;
}

int openmc_next_batch(int* status)
{
  using namespace openmc;
  using openmc::simulation::current_gen;

  // Make sure simulation has been initialized
  if (!simulation::initialized) {
    set_errmsg("Simulation has not been initialized yet.");
    return OPENMC_E_ALLOCATE;
  }

  initialize_batch();


  // =======================================================================
  // LOOP OVER GENERATIONS
  for (current_gen = 1; current_gen <= settings::gen_per_batch; ++current_gen) {

    initialize_generation();

    // Start timer for transport
    simulation::time_transport.start();

    // ====================================================================
    // LOOP OVER PARTICLES
    /*
    #pragma omp parallel for schedule(runtime)
    for (int64_t i_work = 1; i_work <= simulation::work_per_rank; ++i_work) {
      // grab source particle from bank
      Particle p;
      initialize_history(&p, i_work);

      // transport particle
      p.transport();
    }
    */
    /*
    #pragma omp parallel for schedule(runtime)
    for (int64_t i_work = 1; i_work <= simulation::work_per_rank; ++i_work) {
      // grab source particle from bank
      Particle p;
      initialize_history(&p, i_work);

      // transport particle
      p.transport_history_based();
    }
    */

    transport();

    // Accumulate time for transport
    simulation::time_transport.stop();

    finalize_generation();
  }


  finalize_batch();
  

  // Check simulation ending criteria
  if (status) {
    if (simulation::current_batch == settings::n_max_batches) {
      *status = STATUS_EXIT_MAX_BATCH;
    } else if (simulation::satisfy_triggers) {
      *status = STATUS_EXIT_ON_TRIGGER;
    } else {
      *status = STATUS_EXIT_NORMAL;
    }
  }
  return 0;
}

bool openmc_is_statepoint_batch() {
  using namespace openmc;
  using openmc::simulation::current_gen;

  if (!simulation::initialized)
    return false;
  else
    return contains(settings::statepoint_batch, simulation::current_batch);
}

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace simulation {

int current_batch;
int current_gen;
//int64_t current_work;
bool initialized {false};
double keff {1.0};
double keff_std;
double k_col_abs {0.0};
double k_col_tra {0.0};
double k_abs_tra {0.0};
double log_spacing;
int n_lost_particles {0};
bool need_depletion_rx {false};
int restart_batch;
bool satisfy_triggers {false};
int total_gen {0};
double total_weight;
int64_t work_per_rank;

const RegularMesh* entropy_mesh {nullptr};
const RegularMesh* ufs_mesh {nullptr};

std::vector<double> k_generation;
std::vector<int64_t> work_index;

} // namespace simulation

//==============================================================================
// Non-member functions
//==============================================================================

void allocate_banks()
{
  // Allocate source bank
  simulation::source_bank.resize(simulation::work_per_rank);

  if (settings::run_mode == RUN_MODE_EIGENVALUE) {
    simulation::fission_bank.reserve(3*simulation::work_per_rank);
  }
}

void initialize_batch()
{
  // Increment current batch
  ++simulation::current_batch;

  if (settings::run_mode == RUN_MODE_FIXEDSOURCE) {
    int b = simulation::current_batch;
    write_message("Simulating batch " + std::to_string(b), 6);
  }

  // Reset total starting particle weight used for normalizing tallies
  simulation::total_weight = 0.0;

  // Determine if this batch is the first inactive or active batch.
  bool first_inactive = false;
  bool first_active = false;
  if (!settings::restart_run) {
    first_inactive = settings::n_inactive > 0 && simulation::current_batch == 1;
    first_active = simulation::current_batch == settings::n_inactive + 1;
  } else if (simulation::current_batch == simulation::restart_batch + 1){
    first_inactive = simulation::restart_batch < settings::n_inactive;
    first_active = !first_inactive;
  }

  // Manage active/inactive timers and activate tallies if necessary.
  if (first_inactive) {
    simulation::time_inactive.start();
  } else if (first_active) {
    simulation::time_inactive.stop();
    simulation::time_active.start();
    for (auto& t : model::tallies) {
      t->active_ = true;
    }
  }

  // Add user tallies to active tallies list
  setup_active_tallies();
}

void finalize_batch()
{
  // Reduce tallies onto master process and accumulate
  simulation::time_tallies.start();
  accumulate_tallies();
  simulation::time_tallies.stop();

  // Reset global tally results
  if (simulation::current_batch <= settings::n_inactive) {
    xt::view(simulation::global_tallies, xt::all()) = 0.0;
    simulation::n_realizations = 0;
  }

  if (settings::run_mode == RUN_MODE_EIGENVALUE) {
    // Write batch output
    if (mpi::master && settings::verbosity >= 7) print_batch_keff();
  }

  // Check_triggers
  if (mpi::master) check_triggers();
#ifdef OPENMC_MPI
  MPI_Bcast(&simulation::satisfy_triggers, 1, MPI_C_BOOL, 0, mpi::intracomm);
#endif
  if (simulation::satisfy_triggers || (settings::trigger_on &&
      simulation::current_batch == settings::n_max_batches)) {
    settings::statepoint_batch.insert(simulation::current_batch);
  }

  // Write out state point if it's been specified for this batch and is not
  // a CMFD run instance
  if (contains(settings::statepoint_batch, simulation::current_batch)
      && !settings::cmfd_run) {
    if (contains(settings::sourcepoint_batch, simulation::current_batch)
        && settings::source_write && !settings::source_separate) {
      bool b = true;
      openmc_statepoint_write(nullptr, &b);
    } else {
      bool b = false;
      openmc_statepoint_write(nullptr, &b);
    }
  }

  // Write out a separate source point if it's been specified for this batch
  if (contains(settings::sourcepoint_batch, simulation::current_batch)
      && settings::source_write && settings::source_separate) {
    write_source_point(nullptr);
  }

  // Write a continously-overwritten source point if requested.
  if (settings::source_latest) {
    auto filename = settings::path_output + "source.h5";
    write_source_point(filename.c_str());
  }
}

void initialize_generation()
{
  if (settings::run_mode == RUN_MODE_EIGENVALUE) {
    // Clear out the fission bank
    simulation::fission_bank.clear();

    // Count source sites if using uniform fission source weighting
    if (settings::ufs_on) ufs_count_sites();

    // Store current value of tracklength k
    simulation::keff_generation = simulation::global_tallies(
      K_TRACKLENGTH, RESULT_VALUE);
  }
}

/*
struct bank_site_comparator
{
  inline bool operator() (const Particle::Bank & a, const Particle::Bank & b)
  {
    if( a.E < b.E )
      return true;
    else if( a.E > b.E )
      return false;
    else // Energy equal, compare by x-coord
    {
      if(a.r.x < b.r.x )
        return true;
      else if (a.r.x > b.r.x)
        return false;
      else // x-coord equal, compare by y-coord
      {
        if(a.r.x < b.r.x )
          return true;
        else if (a.r.x > b.r.x)
          return false;
        else // y-coord equal, compare by z-coord
        {
          if(a.r.y < b.r.y )
            return true;
          else if (a.r.y > b.r.y)
            return false;
          else // y-coord equal, compare by z-coord
          {
            if(a.r.z < b.r.z )
              return true;
            else if (a.r.z > b.r.z)
              return false;
            else // they are the same
              return false;
          }
        }
      }
    }
  }
};
*/

void finalize_generation()
{
  auto& gt = simulation::global_tallies;

  // Update global tallies with the omp private accumulation variables
  if (settings::run_mode == RUN_MODE_EIGENVALUE) {
    gt(K_COLLISION, RESULT_VALUE) += global_tally_collision;
    gt(K_ABSORPTION, RESULT_VALUE) += global_tally_absorption;
    gt(K_TRACKLENGTH, RESULT_VALUE) += global_tally_tracklength;
  }
  gt(LEAKAGE, RESULT_VALUE) += global_tally_leakage;

  // reset tallies
  if (settings::run_mode == RUN_MODE_EIGENVALUE) {
    global_tally_collision = 0.0;
    global_tally_absorption = 0.0;
    global_tally_tracklength = 0.0;
  }
  global_tally_leakage = 0.0;

  if (settings::run_mode == RUN_MODE_EIGENVALUE) {	
	  // We need to move all the stuff from the shared_fission_bank into the real one.
	  for( int i = 0; i < shared_fission_bank_length; i++ )
		  simulation::fission_bank.push_back(shared_fission_bank[i]);
	  shared_fission_bank_length = 0;

    // Sorts the fission bank so as to allow for reproducibility
    std::stable_sort(simulation::fission_bank.begin(), simulation::fission_bank.end());

    // Distribute fission bank across processors evenly
    synchronize_bank();

    // Calculate shannon entropy
    if (settings::entropy_on) shannon_entropy();

    // Collect results and statistics
    calculate_generation_keff();
    calculate_average_keff();

    // Write generation output
    if (mpi::master && settings::verbosity >= 7) {
      if (simulation::current_gen != settings::gen_per_batch) {
        print_generation();
      }
    }

  } else if (settings::run_mode == RUN_MODE_FIXEDSOURCE) {
    // For fixed-source mode, we need to sample the external source
    fill_source_bank_fixedsource();
  }
}

void initialize_history(Particle* p, int64_t index_source)
{
  // set defaults
  p->from_source(&simulation::source_bank[index_source - 1]);
  p->current_work_ = index_source;

  // Initialize global tally variables
  p->tally_absorption_ =  0.0;
  p->tally_collision_ =   0.0;
  p->tally_tracklength_ = 0.0;
  p->tally_leakage_ =     0.0;

  // set identifier for particle
  p->id_ = simulation::work_index[mpi::rank] + index_source;

  // set random number seed
  int64_t particle_seed = (simulation::total_gen + overall_generation() - 1)
    * settings::n_particles + p->id_;
  init_particle_seeds(particle_seed, p->seeds_);

  // set particle trace
  p->trace_ = false;
  /*
  if (simulation::current_batch == settings::trace_batch &&
      simulation::current_gen == settings::trace_gen &&
      p->id_ == settings::trace_particle) p->trace_ = true;
      */

  // Set particle track.
  p->write_track_ = false;
  if (settings::write_all_tracks) {
    p->write_track_ = true;
  } else if (settings::track_identifiers.size() > 0) {
    for (const auto& t : settings::track_identifiers) {
      if (simulation::current_batch == t[0] &&
          simulation::current_gen == t[1] &&
          p->id_ == t[2]) {
        p->write_track_ = true;
        break;
      }
    }
  }

  // Display message if high verbosity or trace is on
  if (settings::verbosity >= 9 || p->trace_) {
    write_message("Simulating Particle " + std::to_string(p->id_));
  }

  // Add paricle's starting weight to count for normalizing tallies later
  #pragma omp atomic
  simulation::total_weight += p->wgt_;

  // Force calculation of cross-sections by setting last energy to zero
  if (settings::run_CE) {
    for (auto& micro : p->neutron_xs_) micro.last_E = 0.0;
  }

  // Prepare to write out particle track.
  if (p->write_track_) add_particle_track(*p);

  // Every particle starts with no accumulated flux derivative.
  if (!model::active_tallies.empty())
  {
    p->flux_derivs_.resize(model::tally_derivs.size(), 0.0);
    std::fill(p->flux_derivs_.begin(), p->flux_derivs_.end(), 0.0);
  }
  
  // Allocate space for tally filter matches
  p->filter_matches_.resize(model::tally_filters.size());
}

int overall_generation()
{
  using namespace simulation;
  return settings::gen_per_batch*(current_batch - 1) + current_gen;
}

void calculate_work()
{
  // Determine minimum amount of particles to simulate on each processor
  int64_t min_work = settings::n_particles / mpi::n_procs;

  // Determine number of processors that have one extra particle
  int64_t remainder = settings::n_particles % mpi::n_procs;

  int64_t i_bank = 0;
  simulation::work_index.resize(mpi::n_procs + 1);
  simulation::work_index[0] = 0;
  for (int i = 0; i < mpi::n_procs; ++i) {
    // Number of particles for rank i
    int64_t work_i = i < remainder ? min_work + 1 : min_work;

    // Set number of particles
    if (mpi::rank == i) simulation::work_per_rank = work_i;

    // Set index into source bank for rank i
    i_bank += work_i;
    simulation::work_index[i + 1] = i_bank;
  }
}

#ifdef OPENMC_MPI
void broadcast_results() {
  // Broadcast tally results so that each process has access to results
  for (auto& t : model::tallies) {
    // Create a new datatype that consists of all values for a given filter
    // bin and then use that to broadcast. This is done to minimize the
    // chance of the 'count' argument of MPI_BCAST exceeding 2**31
    auto& results = t->results_;

    auto shape = results.shape();
    int count_per_filter = shape[1] * shape[2];
    MPI_Datatype result_block;
    MPI_Type_contiguous(count_per_filter, MPI_DOUBLE, &result_block);
    MPI_Type_commit(&result_block);
    MPI_Bcast(results.data(), shape[0], result_block, 0, mpi::intracomm);
    MPI_Type_free(&result_block);
  }

  // Also broadcast global tally results
  auto& gt = simulation::global_tallies;
  MPI_Bcast(gt.data(), gt.size(), MPI_DOUBLE, 0, mpi::intracomm);

  // These guys are needed so that non-master processes can calculate the
  // combined estimate of k-effective
  double temp[] {simulation::k_col_abs, simulation::k_col_tra,
    simulation::k_abs_tra};
  MPI_Bcast(temp, 3, MPI_DOUBLE, 0, mpi::intracomm);
  simulation::k_col_abs = temp[0];
  simulation::k_col_tra = temp[1];
  simulation::k_abs_tra = temp[2];
}

#endif

void free_memory_simulation()
{
  simulation::k_generation.clear();
  simulation::entropy.clear();
}

} // namespace openmc
