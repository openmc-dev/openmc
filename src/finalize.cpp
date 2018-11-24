#include "openmc/finalize.h"

#include "openmc/capi.h"
#include "openmc/constants.h"
#include "openmc/eigenvalue.h"
#include "openmc/geometry.h"
#include "openmc/message_passing.h"
#include "openmc/nuclide.h"
#include "openmc/random_lcg.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"
#include "openmc/timer.h"
#include "openmc/tallies/tally.h"

using namespace openmc;

// Functions defined in Fortran
extern "C" void free_memory();
extern "C" void reset_timers_f();

int openmc_finalize()
{
  // Clear results
  openmc_reset();

  // Reset global variables
  settings::assume_separate = false;
  settings::check_overlaps = false;
  settings::confidence_intervals = false;
  settings::create_fission_neutrons = true;
  settings::electron_treatment = ELECTRON_LED;
  settings::energy_cutoff = {0.0, 1000.0, 0.0, 0.0};
  settings::entropy_on = false;
  settings::gen_per_batch = 1;
  settings::index_entropy_mesh = -1;
  settings::index_ufs_mesh = -1;
  settings::legendre_to_tabular = true;
  settings::legendre_to_tabular_points = -1;
  settings::n_particles = -1;
  settings::output_summary = true;
  settings::output_tallies = true;
  settings::particle_restart_run = false;
  settings::photon_transport = false;
  settings::reduce_tallies = true;
  settings::res_scat_on = false;
  settings::res_scat_method = RES_SCAT_ARES;
  settings::res_scat_energy_min = 0.01;
  settings::res_scat_energy_max = 1000.0;
  settings::restart_run = false;
  settings::run_CE = true;
  settings::run_mode = -1;
  settings::dagmc = false;
  settings::source_latest = false;
  settings::source_separate = false;
  settings::source_write = true;
  settings::survival_biasing = false;
  settings::temperature_default = 293.6;
  settings::temperature_method = TEMPERATURE_NEAREST;
  settings::temperature_multipole = false;
  settings::temperature_range = {0.0, 0.0};
  settings::temperature_tolerance = 10.0;
  settings::trigger_on = false;
  settings::trigger_predict = false;
  settings::trigger_batch_interval = 1;
  settings::ufs_on = false;
  settings::urr_ptables_on = true;
  settings::verbosity = 7;
  settings::weight_cutoff = 0.25;
  settings::weight_survive = 1.0;
  settings::write_all_tracks = false;
  settings::write_initial_source = false;

  simulation::keff = 1.0;
  simulation::n_lost_particles = 0;
  simulation::satisfy_triggers = false;
  simulation::total_gen = 0;

  data::energy_max = {INFTY, INFTY};
  data::energy_min = {0.0, 0.0};
  n_tallies = 0;
  model::root_universe = -1;
  openmc_set_seed(DEFAULT_SEED);

  // Deallocate arrays
  free_memory();

  // Free all MPI types
#ifdef OPENMC_MPI
  MPI_Type_free(&mpi::bank);
#endif

  return 0;
}

int openmc_reset()
{
  for (int i = 1; i <= n_tallies; ++i) {
    openmc_tally_reset(i);
  }

  // Reset global tallies (can't really use global_tallies() right now because
  // it doesn't have any information about whether the underlying buffer was
  // allocated)
  n_realizations = 0;
  double* buffer = nullptr;
  openmc_global_tallies(&buffer);
  if (buffer) {
    for (int i = 0; i < 3*N_GLOBAL_TALLIES; ++i) {
      buffer[i] = 0.0;
    }
  }

  simulation::k_col_abs = 0.0;
  simulation::k_col_tra = 0.0;
  simulation::k_abs_tra = 0.0;
  simulation::k_sum = {0.0, 0.0};

  // Reset timers
  reset_timers();
  reset_timers_f();

  return 0;
}

int openmc_hard_reset()
{
  // Reset all tallies and timers
  openmc_reset();

  // Reset total generations and keff guess
  simulation::keff = 1.0;
  simulation::total_gen = 0;

  // Reset the random number generator state
  openmc_set_seed(DEFAULT_SEED);
  return 0;
}
