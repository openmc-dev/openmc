#include "openmc/simulation.h"

#include "openmc/bank.h"
#include "openmc/capi.h"
#include "openmc/container_util.h"
#include "openmc/eigenvalue.h"
#include "openmc/error.h"
#include "openmc/event.h"
#include "openmc/geometry_aux.h"
#include "openmc/material.h"
#include "openmc/message_passing.h"
#include "openmc/nuclide.h"
#include "openmc/output.h"
#include "openmc/particle.h"
#include "openmc/photon.h"
#include "openmc/random_lcg.h"
#include "openmc/settings.h"
#include "openmc/source.h"
#include "openmc/state_point.h"
#include "openmc/timer.h"
#include "openmc/tallies/derivative.h"
#include "openmc/tallies/filter.h"
#include "openmc/tallies/tally.h"
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

  // Allocate source and fission banks for eigenvalue simulations
  if (settings::run_mode == RunMode::EIGENVALUE) {
    allocate_banks();
  }

  // If doing an event-based simulation, intialize the particle buffer
  // and event queues
  if (settings::event_based) {
    int64_t event_buffer_length = std::min(simulation::work_per_rank,
      settings::max_particles_in_flight);
    init_event_queues(event_buffer_length);
  }

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
    // Only initialize primary source bank for eigenvalue simulations
    if (settings::run_mode == RunMode::EIGENVALUE) {
      initialize_source();
    }
  }

  // If fixed source and using custom source library then need to load
  if (settings::run_mode == RunMode::FIXED_SOURCE &&
      !settings::path_source_library.empty()) {
    load_custom_source_library();
  }

  // Display header
  if (mpi::master) {
    if (settings::run_mode == RunMode::FIXED_SOURCE) {
      header("FIXED SOURCE TRANSPORT SIMULATION", 3);
    } else if (settings::run_mode == RunMode::EIGENVALUE) {
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

  // If fixed source and using custom source library then need to close
  if (settings::run_mode == RunMode::FIXED_SOURCE &&
      !settings::path_source_library.empty()) {
    close_custom_source_library();
  }

  // Stop timers and show timing statistics
  simulation::time_finalize.stop();
  simulation::time_total.stop();
  if (mpi::master) {
    if (settings::verbosity >= 6) print_runtime();
    if (settings::verbosity >= 4) print_results();
  }
  if (settings::check_overlaps) print_overlap_check();

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

    // Transport loop
    if (settings::event_based) {
      transport_event_based();
    } else {
      transport_history_based();
    }

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
  // Allocate fission bank
  init_fission_bank(3*simulation::work_per_rank);
}

void initialize_batch()
{
  // Increment current batch
  ++simulation::current_batch;

  if (settings::run_mode == RunMode::FIXED_SOURCE) {
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

  if (settings::run_mode == RunMode::EIGENVALUE) {
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
      bool b = (settings::run_mode == RunMode::EIGENVALUE);
      openmc_statepoint_write(nullptr, &b);
    } else {
      bool b = false;
      openmc_statepoint_write(nullptr, &b);
    }
  }

  if (settings::run_mode == RunMode::EIGENVALUE) {
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
}

void initialize_generation()
{
  if (settings::run_mode == RunMode::EIGENVALUE) {
    // Clear out the fission bank
    simulation::fission_bank.resize(0);

    // Count source sites if using uniform fission source weighting
    if (settings::ufs_on) ufs_count_sites();

    // Store current value of tracklength k
    simulation::keff_generation = simulation::global_tallies(
      GlobalTally::K_TRACKLENGTH, TallyResult::VALUE);
  }
}

void finalize_generation()
{
  auto& gt = simulation::global_tallies;

  // Update global tallies with the accumulation variables
  if (settings::run_mode == RunMode::EIGENVALUE) {
    gt(GlobalTally::K_COLLISION, TallyResult::VALUE) += global_tally_collision;
    gt(GlobalTally::K_ABSORPTION, TallyResult::VALUE) += global_tally_absorption;
    gt(GlobalTally::K_TRACKLENGTH, TallyResult::VALUE) += global_tally_tracklength;
  }
  gt(GlobalTally::LEAKAGE, TallyResult::VALUE) += global_tally_leakage;

  // reset tallies
  if (settings::run_mode == RunMode::EIGENVALUE) {
    global_tally_collision = 0.0;
    global_tally_absorption = 0.0;
    global_tally_tracklength = 0.0;
  }
  global_tally_leakage = 0.0;

  if (settings::run_mode == RunMode::EIGENVALUE) {	
    // If using shared memory, stable sort the fission bank (by parent IDs) 
    // so as to allow for reproducibility regardless of which order particles
    // are run in.
    sort_fission_bank();

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

  }
}

void initialize_history(Particle* p, int64_t index_source)
{
  // set defaults
  if (settings::run_mode == RunMode::FIXED_SOURCE) {
    // initialize random number seed
    int64_t id = (simulation::total_gen + overall_generation() - 1)*settings::n_particles +
      simulation::work_index[mpi::rank] + index_source;
    uint64_t seed = init_seed(id, STREAM_SOURCE);
    // sample from external source distribution or custom library then set
    auto site = sample_external_source(&seed);
    p->from_source(&site);
  } else if (settings::run_mode == RunMode::EIGENVALUE) {
    // set defaults for eigenvalue simulations from primary bank
    p->from_source(&simulation::source_bank[index_source - 1]);
  }
  p->current_work_ = index_source;

  // set identifier for particle
  p->id_ = simulation::work_index[mpi::rank] + index_source;

  // set progeny count to zero
  p->n_progeny_ = 0;

  // set random number seed
  int64_t particle_seed = (simulation::total_gen + overall_generation() - 1)
    * settings::n_particles + p->id_;
  init_particle_seeds(particle_seed, p->seeds_);

  // set particle trace
  p->trace_ = false;
  if (simulation::current_batch == settings::trace_batch &&
      simulation::current_gen == settings::trace_gen &&
      p->id_ == settings::trace_particle) p->trace_ = true;

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

  initialize_history_partial(p);
}

void initialize_history_partial(Particle* p)
{
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

void transport_history_based_single_particle(Particle& p)
{
  while (true) {
    p.event_calculate_xs();
    p.event_advance();
    if (p.collision_distance_ > p.boundary_.distance) {
      p.event_cross_surface();
    } else {
      p.event_collide();
    }
    p.event_revive_from_secondary();
    if (!p.alive_)
      break;
  }
  p.event_death();
}

void transport_history_based()
{
  #pragma omp parallel for schedule(runtime)
  for (int64_t i_work = 1; i_work <= simulation::work_per_rank; ++i_work) {
    Particle p;
    initialize_history(&p, i_work);
    transport_history_based_single_particle(p);
  }
}

void transport_event_based()
{
  int64_t remaining_work = simulation::work_per_rank;
  int64_t source_offset = 0;

  // To cap the total amount of memory used to store particle object data, the
  // number of particles in flight at any point in time can bet set. In the case
  // that the maximum in flight particle count is lower than the total number
  // of particles that need to be run this iteration, the event-based transport
  // loop is executed multiple times until all particles have been completed.
  while (remaining_work > 0) {
    // Figure out # of particles to run for this subiteration
    int64_t n_particles = std::min(remaining_work, settings::max_particles_in_flight);

    // Initialize all particle histories for this subiteration
    process_init_events(n_particles, source_offset);

    // Event-based transport loop
    while (true) {
      // Determine which event kernel has the longest queue
      int64_t max = std::max({
        simulation::calculate_fuel_xs_queue.size(),
        simulation::calculate_nonfuel_xs_queue.size(),
        simulation::advance_particle_queue.size(),
        simulation::surface_crossing_queue.size(),
        simulation::collision_queue.size()});

      // Execute event with the longest queue
      if (max == 0) {
        break;
      } else if (max == simulation::calculate_fuel_xs_queue.size()) {
        process_calculate_xs_events(simulation::calculate_fuel_xs_queue);
      } else if (max == simulation::calculate_nonfuel_xs_queue.size()) {
        process_calculate_xs_events(simulation::calculate_nonfuel_xs_queue);
      } else if (max == simulation::advance_particle_queue.size()) {
        process_advance_particle_events();
      } else if (max == simulation::surface_crossing_queue.size()) {
        process_surface_crossing_events();
      } else if (max == simulation::collision_queue.size()) {
        process_collision_events();
      }
    }
    
    // Execute death event for all particles
    process_death_events(n_particles);

    // Adjust remaining work and source offset variables
    remaining_work -= n_particles;
    source_offset += n_particles;
  }
}

} // namespace openmc
