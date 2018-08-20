#include "openmc/particle.h"

#include <algorithm>
#include <sstream>

#include "openmc/capi.h"
#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/hdf5_interface.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"

namespace openmc {

//==============================================================================
// LocalCoord implementation
//==============================================================================

void
LocalCoord::reset()
{
  cell = 0;
  universe = 0;
  lattice = 0;
  lattice_x = 0;
  lattice_y = 0;
  rotated = false;
}

//==============================================================================
// Particle implementation
//==============================================================================

void
Particle::clear()
{
  // reset any coordinate levels
  for (int i=0; i<MAX_COORD; ++i) coord[i].reset();
}

void
Particle::create_secondary(const double* uvw, double E, int type, bool run_CE)
{
  if (n_secondary == MAX_SECONDARY) {
    fatal_error("Too many secondary particles created.");
  }

  int64_t n = n_secondary;
  secondary_bank[n].particle = type;
  secondary_bank[n].wgt = wgt;
  std::copy(coord[0].xyz, coord[0].xyz + 3, secondary_bank[n].xyz);
  std::copy(uvw, uvw + 3, secondary_bank[n].uvw);
  secondary_bank[n].E = E;
  if (!run_CE) secondary_bank[n].E = g;
  n_secondary += 1;
}

void
Particle::initialize()
{
  // Clear coordinate lists
  clear();

  // Set particle to neutron that's alive
  type  = static_cast<int>(ParticleType::neutron);
  alive = true;

  // clear attributes
  surface           = 0;
  cell_born         = 0;
  material          = 0;
  last_material     = 0;
  last_sqrtkT       = 0;
  wgt               = 1.0;
  last_wgt          = 1.0;
  absorb_wgt        = 0.0;
  n_bank            = 0;
  wgt_bank          = 0.0;
  sqrtkT            = -1.0;
  n_collision       = 0;
  fission           = false;
  delayed_group     = 0;
  for (int i=0; i<MAX_DELAYED_GROUPS; ++i) {
    n_delayed_bank[i] = 0;
  }
  g = 0;

  // Set up base level coordinates
  coord[0].universe = C_NONE;
  n_coord = 1;
  last_n_coord = 1;
}

void
Particle::from_source(const Bank* src, bool run_CE, const double* energy_bin_avg)
{
  // set defaults
  initialize();

  // copy attributes from source bank site
  type             = src->particle;
  wgt              = src->wgt;
  last_wgt         = src->wgt;
  std::copy(src->xyz, src->xyz + 3, coord[0].xyz);
  std::copy(src->uvw, src->uvw + 3, coord[0].uvw);
  std::copy(src->xyz, src->xyz + 3, last_xyz_current);
  std::copy(src->xyz, src->xyz + 3, last_xyz);
  std::copy(src->uvw, src->uvw + 3, last_uvw);
  if (run_CE) {
    E = src->E;
    g = 0;
  } else {
    g = static_cast<int>(src->E);
    last_g = static_cast<int>(src->E);
    E = energy_bin_avg[g - 1];
  }
  last_E = E;
}

void
Particle::mark_as_lost(const char* message)
{
  // Print warning and write lost particle file
  warning(message);
  write_restart();

  // Increment number of lost particles
  alive = false;
#pragma omp atomic
  openmc_n_lost_particles += 1;

  // Count the total number of simulated particles (on this processor)
  auto n = openmc_current_batch * gen_per_batch * openmc_work;

  // Abort the simulation if the maximum number of lost particles has been
  // reached
  if (openmc_n_lost_particles >= MAX_LOST_PARTICLES &&
      openmc_n_lost_particles >= REL_MAX_LOST_PARTICLES*n) {
    fatal_error("Maximum number of lost particles has been reached.");
  }
}

void
Particle::write_restart()
{
  // Dont write another restart file if in particle restart mode
  if (openmc_run_mode == RUN_MODE_PARTICLE) return;

  // Set up file name
  std::stringstream filename;
  filename << path_output << "particle_" << openmc_current_batch << '_' << id << ".h5";

#pragma omp critical (WriteParticleRestart)
  {
    // Create file
    hid_t file_id = file_open(filename.str(), 'w');

    // Write filetype and version info
    write_attribute(file_id, "filetype", "particle restart");
    write_attribute(file_id, "version", VERSION_PARTICLE_RESTART);
    write_attribute(file_id, "openmc_version", VERSION);
#ifdef GIT_SHA1
    write_attr_string(file_id, "git_sha1", GIT_SHA1);
#endif

    // Write data to file
    write_dataset(file_id, "current_batch", openmc_current_batch);
    write_dataset(file_id, "generations_per_batch", gen_per_batch);
    write_dataset(file_id, "current_generation", openmc_current_gen);
    write_dataset(file_id, "n_particles", n_particles);
    switch (openmc_run_mode) {
      case RUN_MODE_FIXEDSOURCE:
        write_dataset(file_id, "run_mode", "fixed source");
        break;
      case RUN_MODE_EIGENVALUE:
        write_dataset(file_id, "run_mode", "eigenvalue");
        break;
      case RUN_MODE_PARTICLE:
        write_dataset(file_id, "run_mode", "particle restart");
        break;
    }
    write_dataset(file_id, "id", id);
    write_dataset(file_id, "type", type);

    // Get pointer to source bank
    Bank* src;
    int64_t n;
    openmc_source_bank(&src, &n);

    int64_t i = openmc_current_work;
    write_dataset(file_id, "weight", src[i-1].wgt);
    write_dataset(file_id, "energy", src[i-1].E);
    hsize_t dims[] {3};
    write_double(file_id, 1, dims, "xyz", src[i-1].xyz, false);
    write_double(file_id, 1, dims, "uvw", src[i-1].uvw, false);

    // Close file
    file_close(file_id);
  } // #pragma omp critical
}

//==============================================================================
// Fortran compatibility functions
//==============================================================================

void reset_coord(LocalCoord* c) { c->reset(); }
void particle_clear(Particle* p) { p->clear(); }
void particle_create_secondary(Particle* p, const double* uvw, double E,
                               int type, bool run_CE)
{
  p->create_secondary(uvw, E, type, run_CE);
}
void particle_initialize(Particle* p) { p->initialize(); }
void particle_from_source(Particle* p, const Bank* src, bool run_CE,
                          const double* energy_bin_avg)
{
  p->from_source(src, run_CE, energy_bin_avg);
}
void particle_mark_as_lost(Particle* p, const char* message)
{
  p->mark_as_lost(message);
}
void particle_write_restart(Particle* p) { p->write_restart(); }

} // namespace openmc
