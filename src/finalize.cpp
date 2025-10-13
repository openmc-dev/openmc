#include "openmc/finalize.h"

#include "openmc/bank.h"
#include "openmc/capi.h"
#include "openmc/cmfd_solver.h"
#include "openmc/collision_track.h"
#include "openmc/constants.h"
#include "openmc/cross_sections.h"
#include "openmc/dagmc.h"
#include "openmc/eigenvalue.h"
#include "openmc/event.h"
#include "openmc/geometry.h"
#include "openmc/geometry_aux.h"
#include "openmc/material.h"
#include "openmc/mesh.h"
#include "openmc/message_passing.h"
#include "openmc/nuclide.h"
#include "openmc/photon.h"
#include "openmc/plot.h"
#include "openmc/random_lcg.h"
#include "openmc/random_ray/random_ray_simulation.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"
#include "openmc/source.h"
#include "openmc/surface.h"
#include "openmc/tallies/tally.h"
#include "openmc/thermal.h"
#include "openmc/timer.h"
#include "openmc/volume_calc.h"
#include "openmc/weight_windows.h"

#include "xtensor/xview.hpp"

namespace openmc {

void free_memory()
{
  free_memory_geometry();
  free_memory_surfaces();
  free_memory_material();
  free_memory_volume();
  free_memory_simulation();
  free_memory_photon();
  free_memory_settings();
  free_memory_thermal();
  library_clear();
  nuclides_clear();
  free_memory_source();
  free_memory_mesh();
  free_memory_tally();
  free_memory_bank();
  free_memory_plot();
  free_memory_weight_windows();
  if (mpi::master) {
    free_memory_cmfd();
  }
  if (settings::event_based) {
    free_event_queues();
  }
}

} // namespace openmc

using namespace openmc;

int openmc_finalize()
{
  if (simulation::initialized)
    openmc_simulation_finalize();

  // Clear results
  openmc_reset();

  // Reset timers
  reset_timers();

  // Reset global variables
  settings::assume_separate = false;
  settings::check_overlaps = false;
  settings::confidence_intervals = false;
  settings::create_fission_neutrons = true;
  settings::create_delayed_neutrons = true;
  collision_track::reset_config();
  settings::electron_treatment = ElectronTreatment::LED;
  settings::delayed_photon_scaling = true;
  settings::energy_cutoff = {0.0, 1000.0, 0.0, 0.0};
  settings::time_cutoff = {INFTY, INFTY, INFTY, INFTY};
  settings::entropy_on = false;
  settings::event_based = false;
  settings::free_gas_threshold = 400.0;
  settings::gen_per_batch = 1;
  settings::legendre_to_tabular = true;
  settings::legendre_to_tabular_points = -1;
  settings::material_cell_offsets = true;
  settings::max_lost_particles = 10;
  settings::max_order = 0;
  settings::max_particles_in_flight = 100000;
  settings::max_secondaries = 10000;
  settings::max_particle_events = 1'000'000;
  settings::max_history_splits = 10'000'000;
  settings::max_tracks = 1000;
  settings::max_write_lost_particles = -1;
  settings::n_log_bins = 8000;
  settings::n_inactive = 0;
  settings::n_particles = -1;
  settings::output_summary = true;
  settings::output_tallies = true;
  settings::particle_restart_run = false;
  settings::path_cross_sections.clear();
  settings::path_input.clear();
  settings::path_output.clear();
  settings::path_particle_restart.clear();
  settings::path_sourcepoint.clear();
  settings::path_statepoint.clear();
  settings::photon_transport = false;
  settings::reduce_tallies = true;
  settings::rel_max_lost_particles = 1.0e-6;
  settings::res_scat_on = false;
  settings::res_scat_method = ResScatMethod::rvs;
  settings::res_scat_energy_min = 0.01;
  settings::res_scat_energy_max = 1000.0;
  settings::restart_run = false;
  settings::run_CE = true;
  settings::run_mode = RunMode::UNSET;
  settings::source_latest = false;
  settings::source_rejection_fraction = 0.05;
  settings::source_separate = false;
  settings::source_write = true;
  settings::ssw_cell_id = C_NONE;
  settings::ssw_cell_type = SSWCellType::None;
  settings::ssw_max_particles = 0;
  settings::ssw_max_files = 1;
  settings::survival_biasing = false;
  settings::temperature_default = 293.6;
  settings::temperature_method = TemperatureMethod::NEAREST;
  settings::temperature_multipole = false;
  settings::temperature_range = {0.0, 0.0};
  settings::temperature_tolerance = 10.0;
  settings::trigger_on = false;
  settings::trigger_predict = false;
  settings::trigger_batch_interval = 1;
  settings::uniform_source_sampling = false;
  settings::ufs_on = false;
  settings::urr_ptables_on = true;
  settings::verbosity = 7;
  settings::weight_cutoff = 0.25;
  settings::weight_survive = 1.0;
  settings::weight_windows_file.clear();
  settings::weight_windows_on = false;
  settings::write_all_tracks = false;
  settings::write_initial_source = false;

  simulation::keff = 1.0;
  simulation::need_depletion_rx = false;
  simulation::ssw_current_file = 1;
  simulation::total_gen = 0;

  simulation::entropy_mesh = nullptr;
  simulation::ufs_mesh = nullptr;

  data::energy_max = {INFTY, INFTY};
  data::energy_min = {0.0, 0.0};
  data::temperature_min = 0.0;
  data::temperature_max = INFTY;
  model::root_universe = -1;
  model::plotter_seed = 1;
  openmc::openmc_set_seed(DEFAULT_SEED);
  openmc::openmc_set_stride(DEFAULT_STRIDE);

  // Deallocate arrays
  free_memory();

#ifdef OPENMC_LIBMESH_ENABLED
  settings::libmesh_init.reset();
#endif

  // Free all MPI types
#ifdef OPENMC_MPI
  if (mpi::source_site != MPI_DATATYPE_NULL) {
    MPI_Type_free(&mpi::source_site);
  }
  if (mpi::collision_track_site != MPI_DATATYPE_NULL)
    MPI_Type_free(&mpi::collision_track_site);
#endif

  openmc_reset_random_ray();

  return 0;
}

int openmc_reset()
{

  model::universe_level_counts.clear();

  for (auto& t : model::tallies) {
    t->reset();
  }

  // Reset global tallies
  simulation::n_realizations = 0;
  xt::view(simulation::global_tallies, xt::all()) = 0.0;

  simulation::k_col_abs = 0.0;
  simulation::k_col_tra = 0.0;
  simulation::k_abs_tra = 0.0;
  simulation::k_sum = {0.0, 0.0};
  simulation::satisfy_triggers = false;

  settings::cmfd_run = false;

  simulation::n_lost_particles = 0;

  return 0;
}

int openmc_reset_timers()
{
  reset_timers();
  return 0;
}

int openmc_hard_reset()
{
  // Reset all tallies and timers
  openmc_reset();
  reset_timers();

  // Reset total generations and keff guess
  simulation::keff = 1.0;
  simulation::total_gen = 0;

  // Reset the random number generator state
  openmc::openmc_set_seed(DEFAULT_SEED);
  openmc::openmc_set_stride(DEFAULT_STRIDE);
  return 0;
}
