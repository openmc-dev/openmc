#include "openmc/simulation.h"

#include "openmc/capi.h"
#include "openmc/eigenvalue.h"
#include "openmc/error.h"
#include "openmc/message_passing.h"
#include "openmc/output.h"
#include "openmc/random_lcg.h"
#include "openmc/settings.h"
#include "openmc/source.h"
#include "openmc/timer.h"
#include "openmc/tallies/tally.h"

#include <algorithm>
#include <string>

// data/functions from Fortran side
extern "C" void allocate_banks();
extern "C" void cmfd_init_batch();
extern "C" void cmfd_tally_init();
extern "C" void configure_tallies();
extern "C" void init_tally_routines();
extern "C" void join_bank_from_threads();
extern "C" void load_state_point();
extern "C" void print_columns();
extern "C" void print_generation();
extern "C" void print_results();
extern "C" void print_runtime();
extern "C" void setup_active_tallies();
extern "C" void simulation_init_f();
extern "C" void simulation_finalize_f();
extern "C" void write_tallies();

//==============================================================================
// C API functions
//==============================================================================

// OPENMC_RUN encompasses all the main logic where iterations are performed
// over the batches, generations, and histories in a fixed source or k-eigenvalue
// calculation.

int openmc_run() {
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

  // Set up tally procedure pointers
  init_tally_routines();

  // Allocate source bank, and for eigenvalue simulations also allocate the
  // fission bank
  allocate_banks();

  // Allocate tally results arrays if they're not allocated yet
  configure_tallies();

  // Activate the CMFD tallies
  cmfd_tally_init();

  // Call Fortran initialization
  simulation_init_f();

  // Reset global variables -- this is done before loading state point (as that
  // will potentially populate k_generation and entropy)
  simulation::current_batch = 0;
  simulation::k_generation.clear();
  entropy.clear();
  simulation::need_depletion_rx = false;

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
  time_active.stop();
  time_finalize.start();

  // Deallocate Fortran variables, set tallies to inactive
  simulation_finalize_f();

  // Increment total number of generations
  simulation::total_gen += simulation::current_batch*settings::gen_per_batch;

#ifdef OPENMC_MPI
  broadcast_results()
#endif

  // Write tally results to tallies.out
  if (settings::output_tallies && mpi::master) write_tallies();

  // Deactivate all tallies
  for (int i = 1; i <= n_tallies; ++i) {
    openmc_tally_set_active(i, false);
  }

  // Stop timers and show timing statistics
  time_finalize.stop();
  time_total.stop();
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

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace simulation {

int current_batch;
int current_gen;
int64_t current_work;
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
int64_t work;

std::vector<double> k_generation;
std::vector<int64_t> work_index;

// Threadprivate variables
bool trace;     //!< flag to show debug information
#ifdef _OPENMP
int n_threads {-1};  //!< number of OpenMP threads
int thread_id;  //!< ID of a given thread
#endif

} // namespace simulation

//==============================================================================
// Non-member functions
//==============================================================================

void initialize_batch()
{
  // Increment current batch
  ++simulation::current_batch;

  if (settings::run_mode == RUN_MODE_FIXEDSOURCE) {
    int b = simulation::current_batch;
    write_message("Simulating batch " + std::to_string(b), 6);
  }

  // Reset total starting particle weight used for normalizing tallies
  total_weight = 0.0;

  if ((settings::n_inactive > 0 && simulation::current_batch == 1) ||
      (settings::restart_run && simulation::restart_batch < settings::n_inactive &&
        simulation::current_batch == simulation::restart_batch + 1)) {
    // Turn on inactive timer
    time_inactive.start();
  } else if ((simulation::current_batch == settings::n_inactive + 1) ||
   (settings::restart_run && simulation::restart_batch > settings::n_inactive &&
    simulation::current_batch == simulation::restart_batch + 1)) {
    // Switch from inactive batch timer to active batch timer
    time_inactive.stop();
    time_active.start();

    for (int i = 1; i <= n_tallies; ++i) {
      // TODO: change one-based index
      openmc_tally_set_active(i, true);
    }
  }

  // check CMFD initialize batch
  if (settings::run_mode == RUN_MODE_EIGENVALUE) {
    if (settings::cmfd_run) cmfd_init_batch();
  }

  // Add user tallies to active tallies list
  setup_active_tallies();
}

void initialize_generation()
{
  if (settings::run_mode == RUN_MODE_EIGENVALUE) {
    // Reset number of fission bank sites
    n_bank = 0;

    // Count source sites if using uniform fission source weighting
    if (settings::ufs_on) ufs_count_sites();

    // Store current value of tracklength k
    keff_generation = global_tallies()(K_TRACKLENGTH, RESULT_VALUE);
  }
}

void finalize_generation()
{
  auto gt = global_tallies();

  // Update global tallies with the omp private accumulation variables
#pragma omp parallel
  {
#pragma omp critical(increment_global_tallies)
    {
      if (settings::run_mode == RUN_MODE_EIGENVALUE) {
        gt(K_COLLISION, RESULT_VALUE) += global_tally_collision;
        gt(K_ABSORPTION, RESULT_VALUE) += global_tally_absorption;
        gt(K_TRACKLENGTH, RESULT_VALUE) += global_tally_tracklength;
      }
      gt(LEAKAGE, RESULT_VALUE) += global_tally_leakage;
    }

    // reset threadprivate tallies
    if (settings::run_mode == RUN_MODE_EIGENVALUE) {
      global_tally_collision = 0.0;
      global_tally_absorption = 0.0;
      global_tally_tracklength = 0.0;
    }
    global_tally_leakage = 0.0;
  }


  if (settings::run_mode == RUN_MODE_EIGENVALUE) {
#ifdef _OPENMP
    // Join the fission bank from each thread into one global fission bank
    join_bank_from_threads();
#endif

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
  // Get pointer to source bank
  Bank* source_bank;
  int64_t n;
  openmc_source_bank(&source_bank, &n);

  // set defaults
  p->from_source(&source_bank[index_source - 1]);

  // set identifier for particle
  p->id = simulation::work_index[mpi::rank] + index_source;

  // set random number seed
  int64_t particle_seed = (simulation::total_gen + overall_generation() - 1)
    * settings::n_particles + p->id;
  set_particle_seed(particle_seed);

  // set particle trace
  simulation::trace = false;
  if (simulation::current_batch == settings::trace_batch &&
      simulation::current_gen == settings::trace_gen &&
      p->id == settings::trace_particle) simulation::trace = true;

  // Set particle track.
  p->write_track = false;
  if (settings::write_all_tracks) {
    p->write_track = true;
  } else if (settings::track_identifiers.size() > 0) {
    for (const auto& t : settings::track_identifiers) {
      if (simulation::current_batch == t[0] &&
          simulation::current_gen == t[1] &&
          p->id == t[2]) {
        p->write_track = true;
        break;
      }
    }
  }
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
    if (mpi::rank == i) simulation::work = work_i;

    // Set index into source bank for rank i
    i_bank += work_i;
    simulation::work_index[i + 1] = i_bank;
  }
}

#ifdef OPENMC_MPI
void broadcast_results() {
  // Broadcast tally results so that each process has access to results
  for (int i = 1; i <= n_tallies; ++i) {
    // Create a new datatype that consists of all values for a given filter
    // bin and then use that to broadcast. This is done to minimize the
    // chance of the 'count' argument of MPI_BCAST exceeding 2**31
    auto results = tally_results(i);

    auto shape = results.shape();
    int count_per_filter = shape[1] * shape[2];
    MPI_Datatype result_block;
    MPI_Type_contiguous(count_per_filter, MPI_DOUBLE, &result_block);
    MPI_Type_commit(&result_block);
    MPI_Bcast(results.data(), shape[0], result_block, 0, mpi::intracomm);
    MPI_Type_free(&result_block);
  }

  // Also broadcast global tally results
  auto gt = global_tallies();
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

void broadcast_triggers()
{
  MPI_Bcast(&simulation::satisfy_triggers, 1, MPI_C_BOOL, 0, mpi::intracomm);
}
#endif

//==============================================================================
// Fortran compatibility
//==============================================================================

extern "C" double k_generation(int i) { return simulation::k_generation.at(i - 1); }
extern "C" int k_generation_size() { return simulation::k_generation.size(); }
extern "C" void k_generation_clear() { simulation::k_generation.clear(); }
extern "C" void k_generation_reserve(int i) { simulation::k_generation.reserve(i); }
extern "C" int64_t work_index(int rank) { return simulation::work_index[rank]; }

} // namespace openmc
