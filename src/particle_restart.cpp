#include "openmc/particle_restart.h"

#include "openmc/array.h"
#include "openmc/constants.h"
#include "openmc/hdf5_interface.h"
#include "openmc/mgxs_interface.h"
#include "openmc/nuclide.h"
#include "openmc/output.h"
#include "openmc/particle.h"
#include "openmc/photon.h"
#include "openmc/random_lcg.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"
#include "openmc/tallies/derivative.h"
#include "openmc/tallies/tally.h"
#include "openmc/track_output.h"

#include <algorithm> // for copy
#include <stdexcept>
#include <string>

namespace openmc {

void read_particle_restart(Particle& p, RunMode& previous_run_mode)
{
  // Write meessage
  write_message(
    5, "Loading particle restart file {}", settings::path_particle_restart);

  // Open file
  hid_t file_id = file_open(settings::path_particle_restart, 'r');

  // Read data from file
  read_dataset(file_id, "current_batch", simulation::current_batch);
  read_dataset(file_id, "generations_per_batch", settings::gen_per_batch);
  read_dataset(file_id, "current_generation", simulation::current_gen);
  read_dataset(file_id, "n_particles", settings::n_particles);
  std::string mode;
  read_dataset(file_id, "run_mode", mode);
  if (mode == "eigenvalue") {
    previous_run_mode = RunMode::EIGENVALUE;
  } else if (mode == "fixed source") {
    previous_run_mode = RunMode::FIXED_SOURCE;
  }
  read_dataset(file_id, "id", p.id());
  int type;
  read_dataset(file_id, "type", type);
  p.type() = static_cast<ParticleType>(type);
  read_dataset(file_id, "weight", p.wgt());
  read_dataset(file_id, "energy", p.E());
  read_dataset(file_id, "xyz", p.r());
  read_dataset(file_id, "uvw", p.u());
  read_dataset(file_id, "time", p.time());

  // Set energy group and average energy in multi-group mode
  if (!settings::run_CE) {
    p.g() = p.E();
    p.E() = data::mg.energy_bin_avg_[p.g()];
  }

  // Set particle last attributes
  p.wgt_last() = p.wgt();
  p.r_last_current() = p.r();
  p.r_last() = p.r();
  p.u_last() = p.u();
  p.E_last() = p.E();
  p.g_last() = p.g();
  p.time_last() = p.time();

  // Close hdf5 file
  file_close(file_id);
}

void run_particle_restart()
{
  // Set verbosity high
  settings::verbosity = 10;

  // Initialize nuclear data (energy limits, log grid, etc.)
  initialize_data();

  // Initialize the particle to be tracked
  Particle p;

  // Read in the restart information
  RunMode previous_run_mode;
  read_particle_restart(p, previous_run_mode);

  // write track if that was requested on command line
  if (settings::write_all_tracks) {
    open_track_file();
    p.write_track() = true;
  }

  // Set all tallies to 0 for now (just tracking errors)
  model::tallies.clear();

  // Compute random number seed
  int64_t particle_seed;
  switch (previous_run_mode) {
  case RunMode::EIGENVALUE:
  case RunMode::FIXED_SOURCE:
    particle_seed = (simulation::total_gen + overall_generation() - 1) *
                      settings::n_particles +
                    p.id();
    break;
  default:
    throw std::runtime_error {
      "Unexpected run mode: " +
      std::to_string(static_cast<int>(previous_run_mode))};
  }
  init_particle_seeds(particle_seed, p.seeds());

  // Force calculation of cross-sections by setting last energy to zero
  if (settings::run_CE) {
    p.invalidate_neutron_xs();
  }

  // Prepare to write out particle track.
  if (p.write_track())
    add_particle_track(p);

  // Transport neutron
  transport_history_based_single_particle(p);

  // Write output if particle made it
  print_particle(p);

  if (settings::write_all_tracks) {
    close_track_file();
  }
}

} // namespace openmc
