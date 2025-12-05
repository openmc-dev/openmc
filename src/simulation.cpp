#include "openmc/simulation.h"

#include "openmc/bank.h"
#include "openmc/capi.h"
#include "openmc/collision_track.h"
#include "openmc/container_util.h"
#include "openmc/eigenvalue.h"
#include "openmc/error.h"
#include "openmc/event.h"
#include "openmc/geometry_aux.h"
#include "openmc/ifp.h"
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
#include "openmc/tallies/derivative.h"
#include "openmc/tallies/filter.h"
#include "openmc/tallies/tally.h"
#include "openmc/tallies/trigger.h"
#include "openmc/timer.h"
#include "openmc/track_output.h"
#include "openmc/weight_windows.h"

#ifdef _OPENMP
#include <omp.h>
#endif
#include "xtensor/xview.hpp"

#ifdef OPENMC_MPI
#include <mpi.h>
#endif

#include <fmt/format.h>

#include <algorithm>
#include <cmath>
#include <string>

//==============================================================================
// C API functions
//==============================================================================

// OPENMC_RUN encompasses all the main logic where iterations are performed
// over the batches, generations, and histories in a fixed source or
// k-eigenvalue calculation.

int openmc_run()
{
  openmc::simulation::time_total.start();
  openmc_simulation_init();

  // Ensure that a batch isn't executed in the case that the maximum number of
  // batches has already been run in a restart statepoint file
  int status = 0;
  if (openmc::simulation::current_batch >= openmc::settings::n_max_batches) {
    status = openmc::STATUS_EXIT_MAX_BATCH;
  }

  int err = 0;
  while (status == 0 && err == 0) {
    err = openmc_next_batch(&status);
  }

  openmc_simulation_finalize();
  openmc::simulation::time_total.stop();
  return err;
}

int openmc_simulation_init()
{
  using namespace openmc;

  // Skip if simulation has already been initialized
  if (simulation::initialized)
    return 0;

  // Initialize nuclear data (energy limits, log grid)
  if (settings::run_CE) {
    initialize_data();
  }

  // Determine how much work each process should do
  calculate_work(settings::n_particles);

  // Allocate source, fission and surface source banks.
  allocate_banks();

  // Create track file if needed
  if (!settings::track_identifiers.empty() || settings::write_all_tracks) {
    open_track_file();
  }

  // If doing an event-based simulation, intialize the particle buffer
  // and event queues
  if (settings::event_based) {
    int64_t event_buffer_length =
      std::min(simulation::work_per_rank, settings::max_particles_in_flight);
    init_event_queues(event_buffer_length);
  }

  // Allocate tally results arrays if they're not allocated yet
  for (auto& t : model::tallies) {
    t->set_strides();
    t->init_results();
  }

  // Set up material nuclide index mapping
  for (auto& mat : model::materials) {
    mat->init_nuclide_index();
  }

  // Reset global variables -- this is done before loading state point (as that
  // will potentially populate k_generation and entropy)
  simulation::current_batch = 0;
  simulation::ct_current_file = 1;
  simulation::ssw_current_file = 1;
  simulation::k_generation.clear();
  simulation::entropy.clear();
  openmc_reset();

  // If this is a restart run, load the state point data and binary source
  // file
  if (settings::restart_run) {
    load_state_point();
    write_message("Resuming simulation...", 6);
  } else {
    // Only initialize primary source bank for eigenvalue simulations
    if (settings::run_mode == RunMode::EIGENVALUE &&
        settings::solver_type == SolverType::MONTE_CARLO) {
      initialize_source();
    }
  }

  // Display header
  if (mpi::master) {
    if (settings::run_mode == RunMode::FIXED_SOURCE) {
      if (settings::solver_type == SolverType::MONTE_CARLO) {
        header("FIXED SOURCE TRANSPORT SIMULATION", 3);
      } else if (settings::solver_type == SolverType::RANDOM_RAY) {
        header("FIXED SOURCE TRANSPORT SIMULATION (RANDOM RAY SOLVER)", 3);
      }
    } else if (settings::run_mode == RunMode::EIGENVALUE) {
      if (settings::solver_type == SolverType::MONTE_CARLO) {
        header("K EIGENVALUE SIMULATION", 3);
      } else if (settings::solver_type == SolverType::RANDOM_RAY) {
        header("K EIGENVALUE SIMULATION (RANDOM RAY SOLVER)", 3);
      }
      if (settings::verbosity >= 7)
        print_columns();
    }
  }

  // load weight windows from file
  if (!settings::weight_windows_file.empty()) {
    openmc_weight_windows_import(settings::weight_windows_file.c_str());
  }

  // Set flag indicating initialization is done
  simulation::initialized = true;
  return 0;
}

int openmc_simulation_finalize()
{
  using namespace openmc;

  // Skip if simulation was never run
  if (!simulation::initialized)
    return 0;

  // Stop active batch timer and start finalization timer
  simulation::time_active.stop();
  simulation::time_finalize.start();

  // Clear material nuclide mapping
  for (auto& mat : model::materials) {
    mat->mat_nuclide_index_.clear();
  }

  // Close track file if open
  if (!settings::track_identifiers.empty() || settings::write_all_tracks) {
    close_track_file();
  }

  // Increment total number of generations
  simulation::total_gen += simulation::current_batch * settings::gen_per_batch;

#ifdef OPENMC_MPI
  broadcast_results();
#endif

  // Write tally results to tallies.out
  if (settings::output_tallies && mpi::master)
    write_tallies();

  // If weight window generators are present in this simulation,
  // write a weight windows file
  if (variance_reduction::weight_windows_generators.size() > 0) {
    openmc_weight_windows_export();
  }

  // Deactivate all tallies
  for (auto& t : model::tallies) {
    t->active_ = false;
  }

  // Stop timers and show timing statistics
  simulation::time_finalize.stop();
  simulation::time_total.stop();
  if (mpi::master) {
    if (settings::solver_type != SolverType::RANDOM_RAY) {
      if (settings::verbosity >= 6)
        print_runtime();
      if (settings::verbosity >= 4)
        print_results();
    }
  }
  if (settings::check_overlaps)
    print_overlap_check();

  // Reset flags
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
      if (settings::use_shared_secondary_bank) {
        transport_history_based_shared_secondary();
      } else {
        transport_history_based();
      }
    }

    // Accumulate time for transport
    simulation::time_transport.stop();

    finalize_generation();
  }

  finalize_batch();

  // Check simulation ending criteria
  if (status) {
    if (simulation::current_batch >= settings::n_max_batches) {
      *status = STATUS_EXIT_MAX_BATCH;
    } else if (simulation::satisfy_triggers) {
      *status = STATUS_EXIT_ON_TRIGGER;
    } else {
      *status = STATUS_EXIT_NORMAL;
    }
  }
  return 0;
}

bool openmc_is_statepoint_batch()
{
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

int ct_current_file;
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
int ssw_current_file;
int total_gen {0};
double total_weight;
int64_t work_per_rank;

const RegularMesh* entropy_mesh {nullptr};
const RegularMesh* ufs_mesh {nullptr};

vector<double> k_generation;
vector<int64_t> work_index;

int64_t simulation_particles_completed {0};

} // namespace simulation

//==============================================================================
// Non-member functions
//==============================================================================

void allocate_banks()
{
  if (settings::run_mode == RunMode::EIGENVALUE &&
      settings::solver_type == SolverType::MONTE_CARLO) {
    // Allocate source bank
    simulation::source_bank.resize(simulation::work_per_rank);

    // Allocate fission bank
    init_fission_bank(3 * simulation::work_per_rank);

    // Allocate IFP bank
    if (settings::ifp_on) {
      resize_simulation_ifp_banks();
    }
  }

  if (settings::surf_source_write) {
    // Allocate surface source bank
    simulation::surf_source_bank.reserve(settings::ssw_max_particles);
  }

  if (settings::collision_track) {
    // Allocate collision track bank
    collision_track_reserve_bank();
  }
}

void initialize_batch()
{
  // Increment current batch
  ++simulation::current_batch;
  if (settings::run_mode == RunMode::FIXED_SOURCE) {
    if (settings::solver_type == SolverType::RANDOM_RAY &&
        simulation::current_batch < settings::n_inactive + 1) {
      write_message(
        6, "Simulating batch {:<4} (inactive)", simulation::current_batch);
    } else {
      write_message(6, "Simulating batch {}", simulation::current_batch);
    }
  }

  // Reset total starting particle weight used for normalizing tallies
  simulation::total_weight = 0.0;

  // Determine if this batch is the first inactive or active batch.
  bool first_inactive = false;
  bool first_active = false;
  if (!settings::restart_run) {
    first_inactive = settings::n_inactive > 0 && simulation::current_batch == 1;
    first_active = simulation::current_batch == settings::n_inactive + 1;
  } else if (simulation::current_batch == simulation::restart_batch + 1) {
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

  // update weight windows if needed
  for (const auto& wwg : variance_reduction::weight_windows_generators) {
    wwg->update();
  }

  // Reset global tally results
  if (simulation::current_batch <= settings::n_inactive) {
    xt::view(simulation::global_tallies, xt::all()) = 0.0;
    simulation::n_realizations = 0;
  }

  // Check_triggers
  if (mpi::master)
    check_triggers();
#ifdef OPENMC_MPI
  MPI_Bcast(&simulation::satisfy_triggers, 1, MPI_C_BOOL, 0, mpi::intracomm);
#endif
  if (simulation::satisfy_triggers ||
      (settings::trigger_on &&
        simulation::current_batch == settings::n_max_batches)) {
    settings::statepoint_batch.insert(simulation::current_batch);
  }

  // Write out state point if it's been specified for this batch and is not
  // a CMFD run instance
  if (contains(settings::statepoint_batch, simulation::current_batch) &&
      !settings::cmfd_run) {
    if (contains(settings::sourcepoint_batch, simulation::current_batch) &&
        settings::source_write && !settings::source_separate) {
      bool b = (settings::run_mode == RunMode::EIGENVALUE);
      openmc_statepoint_write(nullptr, &b);
    } else {
      bool b = false;
      openmc_statepoint_write(nullptr, &b);
    }
  }

  if (settings::run_mode == RunMode::EIGENVALUE) {
    // Write out a separate source point if it's been specified for this batch
    if (contains(settings::sourcepoint_batch, simulation::current_batch) &&
        settings::source_write && settings::source_separate) {

      // Determine width for zero padding
      int w = std::to_string(settings::n_max_batches).size();
      std::string source_point_filename = fmt::format("{0}source.{1:0{2}}",
        settings::path_output, simulation::current_batch, w);
      span<SourceSite> bankspan(simulation::source_bank);
      write_source_point(source_point_filename, bankspan,
        simulation::work_index, settings::source_mcpl_write);
    }

    // Write a continously-overwritten source point if requested.
    if (settings::source_latest) {
      auto filename = settings::path_output + "source";
      span<SourceSite> bankspan(simulation::source_bank);
      write_source_point(filename, bankspan, simulation::work_index,
        settings::source_mcpl_write);
    }
  }

  // Write out surface source if requested.
  if (settings::surf_source_write &&
      simulation::ssw_current_file <= settings::ssw_max_files) {
    bool last_batch = (simulation::current_batch == settings::n_batches);
    if (simulation::surf_source_bank.full() || last_batch) {
      // Determine appropriate filename
      auto filename = fmt::format("{}surface_source.{}", settings::path_output,
        simulation::current_batch);
      if (settings::ssw_max_files == 1 ||
          (simulation::ssw_current_file == 1 && last_batch)) {
        filename = settings::path_output + "surface_source";
      }

      // Get span of source bank and calculate parallel index vector
      auto surf_work_index = mpi::calculate_parallel_index_vector(
        simulation::surf_source_bank.size());
      span<SourceSite> surfbankspan(simulation::surf_source_bank.begin(),
        simulation::surf_source_bank.size());

      // Write surface source file
      write_source_point(
        filename, surfbankspan, surf_work_index, settings::surf_mcpl_write);

      // Reset surface source bank and increment counter
      simulation::surf_source_bank.clear();
      if (!last_batch && settings::ssw_max_files >= 1) {
        simulation::surf_source_bank.reserve(settings::ssw_max_particles);
      }
      ++simulation::ssw_current_file;
    }
  }
  // Write collision track file if requested
  if (settings::collision_track) {
    collision_track_flush_bank();
  }
}

void initialize_generation()
{
  if (settings::run_mode == RunMode::EIGENVALUE) {
    // Clear out the fission bank
    simulation::fission_bank.resize(0);

    // Count source sites if using uniform fission source weighting
    if (settings::ufs_on)
      ufs_count_sites();

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
    gt(GlobalTally::K_ABSORPTION, TallyResult::VALUE) +=
      global_tally_absorption;
    gt(GlobalTally::K_TRACKLENGTH, TallyResult::VALUE) +=
      global_tally_tracklength;
  }
  gt(GlobalTally::LEAKAGE, TallyResult::VALUE) += global_tally_leakage;

  // reset tallies
  if (settings::run_mode == RunMode::EIGENVALUE) {
    global_tally_collision = 0.0;
    global_tally_absorption = 0.0;
    global_tally_tracklength = 0.0;
  }
  global_tally_leakage = 0.0;

  if (settings::run_mode == RunMode::EIGENVALUE &&
      settings::solver_type == SolverType::MONTE_CARLO) {
    // If using shared memory, stable sort the fission bank (by parent IDs)
    // so as to allow for reproducibility regardless of which order particles
    // are run in.
    sort_bank(simulation::fission_bank, true);

    // Distribute fission bank across processors evenly
    synchronize_bank();
  }

  if (settings::run_mode == RunMode::EIGENVALUE) {

    // Calculate shannon entropy
    if (settings::entropy_on &&
        settings::solver_type == SolverType::MONTE_CARLO)
      shannon_entropy();

    // Collect results and statistics
    calculate_generation_keff();
    calculate_average_keff();

    // Write generation output
    if (mpi::master && settings::verbosity >= 7) {
      print_generation();
    }
  }
}

void sample_particle(Particle& p, int64_t index_source)
{
  // Sample a particle from the source bank
  if (settings::run_mode == RunMode::EIGENVALUE) {
    p.from_source(&simulation::source_bank[index_source - 1]);
  } else if (settings::run_mode == RunMode::FIXED_SOURCE) {
    // initialize random number seed
    int64_t id;
    if (settings::use_shared_secondary_bank) {
      id = simulation::work_index[mpi::rank] + index_source +
           simulation::simulation_particles_completed;
    } else {
      id = (simulation::total_gen + overall_generation() - 1) *
             settings::n_particles +
           simulation::work_index[mpi::rank] + index_source;
    }
    uint64_t seed = init_seed(id, STREAM_SOURCE);
    // sample from external source distribution or custom library then set
    auto site = sample_external_source(&seed);
    p.from_source(&site);
  }
}

void initialize_history(Particle& p, int64_t index_source, bool is_secondary)
{
  // set defaults
  if (!is_secondary) {
    sample_particle(p, index_source);
  }

  p.current_work() = index_source - 1;

  // set identifier for particle
  if (settings::use_shared_secondary_bank) {
    p.id() = simulation::work_index[mpi::rank] + index_source +
             simulation::simulation_particles_completed;
  } else {
    p.id() = simulation::work_index[mpi::rank] + index_source;
  }

  // set progeny count to zero
  p.n_progeny() = 0;

  // Reset particle event counter
  p.n_event() = 0;

  // Reset split counter
  p.n_split() = 0;

  // Reset weight window ratio
  p.ww_factor() = 0.0;

  // set particle history start weight
  p.wgt_born() = p.wgt();

  // Reset pulse_height_storage
  std::fill(p.pht_storage().begin(), p.pht_storage().end(), 0);

  // set random number seed
  int64_t particle_seed;
  if (settings::use_shared_secondary_bank) {
    particle_seed = p.id() - 1;
  } else {
    particle_seed = (simulation::total_gen + overall_generation() - 1) *
                      settings::n_particles +
                    p.id();
  }
  init_particle_seeds(particle_seed, p.seeds());

  // set particle trace
  p.trace() = false;
  if (simulation::current_batch == settings::trace_batch &&
      simulation::current_gen == settings::trace_gen &&
      p.id() == settings::trace_particle)
    p.trace() = true;

  // Set particle track.
  p.write_track() = check_track_criteria(p);

  // Set the particle's initial weight window value.
  if (!is_secondary) {
    p.wgt_ww_born() = -1.0;
    apply_weight_windows(p);
  }

  // Display message if high verbosity or trace is on
  if (settings::verbosity >= 9 || p.trace()) {
    write_message("Simulating Particle {}", p.id());
  }

  // Add particle's starting weight to count for normalizing tallies later
  if (!is_secondary) {
#pragma omp atomic
    simulation::total_weight += p.wgt();
  }

  // Force calculation of cross-sections by setting last energy to zero
  if (settings::run_CE) {
    p.invalidate_neutron_xs();
  }

  // Prepare to write out particle track.
  if (p.write_track())
    add_particle_track(p);
}

int overall_generation()
{
  using namespace simulation;
  return settings::gen_per_batch * (current_batch - 1) + current_gen;
}

void calculate_work(int64_t n_particles)
{
  // Determine minimum amount of particles to simulate on each processor
  int64_t min_work = n_particles / mpi::n_procs;

  // Determine number of processors that have one extra particle
  int64_t remainder = n_particles % mpi::n_procs;

  int64_t i_bank = 0;
  simulation::work_index.resize(mpi::n_procs + 1);
  simulation::work_index[0] = 0;
  for (int i = 0; i < mpi::n_procs; ++i) {
    // Number of particles for rank i
    int64_t work_i = i < remainder ? min_work + 1 : min_work;

    // Set number of particles
    if (mpi::rank == i)
      simulation::work_per_rank = work_i;

    // Set index into source bank for rank i
    i_bank += work_i;
    simulation::work_index[i + 1] = i_bank;
  }
}

void initialize_data()
{
  // Determine minimum/maximum energy for incident neutron/photon data
  data::energy_max = {INFTY, INFTY, INFTY, INFTY};
  data::energy_min = {0.0, 0.0, 0.0, 0.0};

  for (const auto& nuc : data::nuclides) {
    if (nuc->grid_.size() >= 1) {
      int neutron = static_cast<int>(ParticleType::neutron);
      data::energy_min[neutron] =
        std::max(data::energy_min[neutron], nuc->grid_[0].energy.front());
      data::energy_max[neutron] =
        std::min(data::energy_max[neutron], nuc->grid_[0].energy.back());
    }
  }

  if (settings::photon_transport) {
    for (const auto& elem : data::elements) {
      if (elem->energy_.size() >= 1) {
        int photon = static_cast<int>(ParticleType::photon);
        int n = elem->energy_.size();
        data::energy_min[photon] =
          std::max(data::energy_min[photon], std::exp(elem->energy_(1)));
        data::energy_max[photon] =
          std::min(data::energy_max[photon], std::exp(elem->energy_(n - 1)));
      }
    }

    if (settings::electron_treatment == ElectronTreatment::TTB) {
      // Determine if minimum/maximum energy for bremsstrahlung is greater/less
      // than the current minimum/maximum
      if (data::ttb_e_grid.size() >= 1) {
        int photon = static_cast<int>(ParticleType::photon);
        int electron = static_cast<int>(ParticleType::electron);
        int positron = static_cast<int>(ParticleType::positron);
        int n_e = data::ttb_e_grid.size();

        const std::vector<int> charged = {electron, positron};
        for (auto t : charged) {
          data::energy_min[t] = std::exp(data::ttb_e_grid(1));
          data::energy_max[t] = std::exp(data::ttb_e_grid(n_e - 1));
        }

        data::energy_min[photon] =
          std::max(data::energy_min[photon], data::energy_min[electron]);

        data::energy_max[photon] =
          std::min(data::energy_max[photon], data::energy_max[electron]);
      }
    }
  }

  // Show which nuclide results in lowest energy for neutron transport
  for (const auto& nuc : data::nuclides) {
    // If a nuclide is present in a material that's not used in the model, its
    // grid has not been allocated
    if (nuc->grid_.size() > 0) {
      double max_E = nuc->grid_[0].energy.back();
      int neutron = static_cast<int>(ParticleType::neutron);
      if (max_E == data::energy_max[neutron]) {
        write_message(7, "Maximum neutron transport energy: {} eV for {}",
          data::energy_max[neutron], nuc->name_);
        if (mpi::master && data::energy_max[neutron] < 20.0e6) {
          warning("Maximum neutron energy is below 20 MeV. This may bias "
                  "the results.");
        }
        break;
      }
    }
  }

  // Set up logarithmic grid for nuclides
  for (auto& nuc : data::nuclides) {
    nuc->init_grid();
  }
  int neutron = static_cast<int>(ParticleType::neutron);
  simulation::log_spacing =
    std::log(data::energy_max[neutron] / data::energy_min[neutron]) /
    settings::n_log_bins;
}

#ifdef OPENMC_MPI
void broadcast_results()
{
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
  double temp[] {
    simulation::k_col_abs, simulation::k_col_tra, simulation::k_abs_tra};
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
  while (p.alive()) {
    p.event_calculate_xs();
    if (p.alive()) {
      p.event_advance();
    }
    if (p.alive()) {
      if (p.collision_distance() > p.boundary().distance()) {
        p.event_cross_surface();
      } else if (p.alive()) {
        p.event_collide();
      }
    }
    // If particle has too many events, display warning and kill it
    p.n_event()++;
    if (p.n_event() == settings::max_particle_events) {
      warning("Particle " + std::to_string(p.id()) +
              " underwent maximum number of events.");
      p.wgt() = 0.0;
    }
    if (!p.alive() && !settings::use_shared_secondary_bank &&
        !p.local_secondary_bank().empty()) {
      SourceSite& site = p.local_secondary_bank().back();
      p.event_revive_from_secondary(site);
      p.local_secondary_bank().pop_back();
    }
  }
  p.event_death();
}

void transport_history_based()
{
#pragma omp parallel for schedule(runtime)
  for (int64_t i_work = 1; i_work <= simulation::work_per_rank; ++i_work) {
    Particle p;
    initialize_history(p, i_work, false);
    transport_history_based_single_particle(p);
  }
}

void debug_validate_local_bank_ordering(
  const SharedArray<SourceSite>& local_bank)
{
  for (int64_t i = 1; i < local_bank.size(); i++) {
    const SourceSite& a = local_bank[i - 1];
    const SourceSite& b = local_bank[i];
    if (a.parent_id > b.parent_id ||
        (a.parent_id == b.parent_id && a.progeny_id > b.progeny_id)) {
          fmt::print("Error at local_bank[{}]: ({}, {}) > ({}, {})\n", i,
            a.parent_id, a.progeny_id, b.parent_id, b.progeny_id);
      fatal_error("Local secondary bank is not properly ordered by parent and "
                  "progeny IDs.");
    }
  }
}

// Another helpfer function that validates if the bank is GLOBALLY sorted across
// all MPI ranks. This is an expensive check and should only be used for
// debugging purposes.
void debug_validate_global_bank_ordering(
  const SharedArray<SourceSite>& local_bank)
{
#ifdef OPENMC_MPI
  // Each rank (except rank 0) receives the last entry from the previous rank
  // and compares it to its first entry to ensure global ordering.

  // Data to send/receive: parent_id and progeny_id
  int64_t last_entry[2] = {0, 0};  // parent_id, progeny_id
  int64_t recv_entry[2] = {0, 0};

  // If this rank has entries, get the last one
  if (local_bank.size() > 0) {
    const SourceSite& last = local_bank[local_bank.size() - 1];
    last_entry[0] = last.parent_id;
    last_entry[1] = last.progeny_id;
  }

  // Send last entry to next rank, receive from previous rank
  MPI_Request send_request, recv_request;
  MPI_Status status;

  // Rank n sends its last entry to rank n+1
  if (mpi::rank < mpi::n_procs - 1) {
    MPI_Isend(last_entry, 2, MPI_INT64_T, mpi::rank + 1, 0, mpi::intracomm,
      &send_request);
  }

  // Rank n receives from rank n-1
  if (mpi::rank > 0) {
    MPI_Irecv(recv_entry, 2, MPI_INT64_T, mpi::rank - 1, 0, mpi::intracomm,
      &recv_request);
  }

  // Wait for communication to complete
  if (mpi::rank < mpi::n_procs - 1) {
    MPI_Wait(&send_request, &status);
  }
  if (mpi::rank > 0) {
    MPI_Wait(&recv_request, &status);
  }

  // Now validate: if this rank has entries and received data from previous rank,
  // check that our first entry is >= the last entry from previous rank
  if (mpi::rank > 0 && local_bank.size() > 0) {
    const SourceSite& first = local_bank[0];
    int64_t prev_parent_id = recv_entry[0];
    int64_t prev_progeny_id = recv_entry[1];

    // Check ordering: previous rank's last entry should be <= this rank's first
    if (prev_parent_id > first.parent_id ||
        (prev_parent_id == first.parent_id &&
         prev_progeny_id > first.progeny_id)) {
      fatal_error(fmt::format(
        "Global secondary bank ordering violated between rank {} and {}: "
        "prev_last=({}, {}), curr_first=({}, {})",
        mpi::rank - 1, mpi::rank, prev_parent_id, prev_progeny_id,
        first.parent_id, first.progeny_id));
    }
  }

  // Synchronize all ranks before continuing
  MPI_Barrier(mpi::intracomm);
#endif
}

// The shared secondary bank transport algorithm works in two phases. In the
// first phase, all primary particles are sampled then transported, and their
// secondary particles are deposited into a shared secondary bank. The second
// phase occurs in a loop, where all secondary particles in the shared secondary
// bank are transported. Any secondary particles generated during this phase are
// deposited back into the shared secondary bank. The shared secondary bank is
// sorted for consistent ordering and load balanced across MPI ranks. This loop
// continues until there are no more secondary particles left to transport.
void transport_history_based_shared_secondary()
{
  // Shared secondary banks for reading and writing
  SharedArray<SourceSite> shared_secondary_bank_read;
  SharedArray<SourceSite> shared_secondary_bank_write;

  if (mpi::master) {
    write_message(fmt::format(" Primogenitor            particles: {}",
                    settings::n_particles),
      6);
  }

  simulation::progeny_per_particle.resize(simulation::work_per_rank);

  // Phase 1: Transport primary particles and deposit first generation of
  // secondaries in the shared secondary bank
#pragma omp parallel for schedule(runtime)
  for (int64_t i = 1; i <= simulation::work_per_rank; i++) {
    Particle p;
    initialize_history(p, i, false);
    transport_history_based_single_particle(p);

    // Transfer all secondary particles to the shared secondary bank
#pragma omp critical(shared_secondary_bank)
    {
      for (auto& site : p.local_secondary_bank()) {
        shared_secondary_bank_write.thread_unsafe_append(site);
      }
    }
  }

  simulation::simulation_particles_completed += settings::n_particles;

  // Phase 2: Now that the secondary bank has been populated, enter loop over
  // all secondary generations
  int n_generation_depth = 1;
  int64_t alive_secondary = 1;
  while (alive_secondary) {

    // Naive serial O(nlogn) sort of the shared secondary bank:
    // Order the shared secondary bank by parent ID then progeny
    // ID to ensure reproducibility.
    /*
    std::sort(shared_secondary_bank_write.begin(),
      shared_secondary_bank_write.end(),
      [](const SourceSite& a, const SourceSite& b) {
        if (a.parent_id != b.parent_id) {
          return a.parent_id < b.parent_id;
        } else {
          return a.progeny_id < b.progeny_id;
        }
      });
    */
    sort_bank(shared_secondary_bank_write, false);

    // Debugging: Validate that the bank is sorted:
    debug_validate_local_bank_ordering(shared_secondary_bank_write);

    // Synchronize the shared secondary bank amongst all MPI ranks, such
    // that each MPI rank has an approximately equal number of secondary
    // particles. Also reports the total number of secondaries alive across
    // all MPI ranks.
    alive_secondary =
      synchronize_global_secondary_bank(shared_secondary_bank_write);
    
    debug_validate_local_bank_ordering(shared_secondary_bank_write);
    debug_validate_global_bank_ordering(shared_secondary_bank_write);

    // Recalculate work for each MPI rank based on number of alive secondary
    // particles
    calculate_work(alive_secondary);

    // Display the number of secondary particles in this generation. This
    // is useful for user monitoring so as to see if the secondary population is
    // exploding and to determine how many generations of secondaries are being
    // transported.
    if (mpi::master) {
      write_message(fmt::format(" Secondary generation {:<2} particles: {}",
                      n_generation_depth, alive_secondary),
        6);
    }

    shared_secondary_bank_read = std::move(shared_secondary_bank_write);
    shared_secondary_bank_write = SharedArray<SourceSite>();
    simulation::progeny_per_particle.resize(shared_secondary_bank_read.size());

    // Transport all secondary particles from the shared secondary bank
#pragma omp parallel for schedule(runtime)
    for (int64_t i = 1; i <= shared_secondary_bank_read.size(); i++) {
      Particle p;
      initialize_history(p, i, true);
      SourceSite& site = shared_secondary_bank_read[i-1];
      p.event_revive_from_secondary(site);
      if (p.alive()) {
        transport_history_based_single_particle(p);
      }
#pragma omp critical(shared_secondary_bank)
      {
        for (auto& site : p.local_secondary_bank()) {
          shared_secondary_bank_write.thread_unsafe_append(site);
        }
      }
    } // End of transport loop over particles in shared secondary bank
    n_generation_depth++;
    simulation::simulation_particles_completed += alive_secondary;
  } // End of loop over secondary generations
  
  // Reset work so that fission bank etc works correctly
  calculate_work(settings::n_particles);
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
    int64_t n_particles =
      std::min(remaining_work, settings::max_particles_in_flight);

    // Initialize all particle histories for this subiteration
    process_init_events(n_particles, source_offset);

    // Event-based transport loop
    while (true) {
      // Determine which event kernel has the longest queue
      int64_t max = std::max({simulation::calculate_fuel_xs_queue.size(),
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
