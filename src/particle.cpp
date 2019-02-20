#include "openmc/particle.h"

#include <algorithm> // copy, min
#include <cmath>     // log, abs, copysign
#include <sstream>

#include "openmc/bank.h"
#include "openmc/capi.h"
#include "openmc/cell.h"
#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/geometry.h"
#include "openmc/hdf5_interface.h"
#include "openmc/material.h"
#include "openmc/message_passing.h"
#include "openmc/mgxs_interface.h"
#include "openmc/nuclide.h"
#include "openmc/physics.h"
#include "openmc/physics_mg.h"
#include "openmc/random_lcg.h"
#include "openmc/settings.h"
#include "openmc/surface.h"
#include "openmc/simulation.h"
#include "openmc/tallies/derivative.h"
#include "openmc/tallies/tally.h"
#include "openmc/tallies/tally_scoring.h"
#include "openmc/track_output.h"

namespace openmc {

//==============================================================================
// LocalCoord implementation
//==============================================================================

void
LocalCoord::reset()
{
  cell = C_NONE;
  universe = C_NONE;
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
  cell_born         = C_NONE;
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
Particle::from_source(const Bank* src)
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
  if (settings::run_CE) {
    E = src->E;
    g = 0;
  } else {
    g = static_cast<int>(src->E);
    last_g = static_cast<int>(src->E);
    E = data::energy_bin_avg[g - 1];
  }
  last_E = E;
}

void
Particle::transport()
{
  // Display message if high verbosity or trace is on
  if (settings::verbosity >= 9 || simulation::trace) {
     write_message("Simulating Particle " + std::to_string(id));
  }

  // Initialize number of events to zero
  int n_event = 0;

  // Add paricle's starting weight to count for normalizing tallies later
  #pragma omp atomic
  simulation::total_weight += wgt;

  // Force calculation of cross-sections by setting last energy to zero
  if (settings::run_CE) {
    for (int i = 0; i < data::nuclides.size(); ++i) {
      simulation::micro_xs[i].last_E = 0.0;
    }
  }

  // Prepare to write out particle track.
  if (write_track) add_particle_track();

  // Every particle starts with no accumulated flux derivative.
  if (!model::active_tallies.empty()) zero_flux_derivs();

  while (true) {
    // Set the random number stream
    if (type == static_cast<int>(ParticleType::neutron)) {
      prn_set_stream(STREAM_TRACKING);
    } else {
      prn_set_stream(STREAM_PHOTON);
    }

    // Store pre-collision particle properties
    last_wgt = wgt;
    last_E = E;
    std::copy(coord[0].uvw, coord[0].uvw + 3, last_uvw);
    std::copy(coord[0].xyz, coord[0].xyz + 3, last_xyz);

    // If the cell hasn't been determined based on the particle's location,
    // initiate a search for the current cell. This generally happens at the
    // beginning of the history and again for any secondary particles
    if (coord[n_coord - 1].cell == C_NONE) {
      if (!find_cell(this, false)) {
        this->mark_as_lost("Could not find the cell containing particle "
          + std::to_string(id));
        return;
      }

      // set birth cell attribute
      if (cell_born == C_NONE) cell_born = coord[n_coord - 1].cell;
    }

    // Write particle track.
    if (write_track) write_particle_track(*this);

    if (settings::check_overlaps) check_cell_overlap(this);

    // Calculate microscopic and macroscopic cross sections
    if (material != MATERIAL_VOID) {
      if (settings::run_CE) {
        if (material != last_material || sqrtkT != last_sqrtkT) {
          // If the material is the same as the last material and the
          // temperature hasn't changed, we don't need to lookup cross
          // sections again.
          model::materials[material - 1]->calculate_xs(*this);
        }
      } else {
        // Get the MG data
        calculate_xs_c(material, g, sqrtkT, coord[n_coord-1].uvw,
          simulation::material_xs.total, simulation::material_xs.absorption,
          simulation::material_xs.nu_fission);

        // Finally, update the particle group while we have already checked
        // for if multi-group
        last_g = g;
      }
    } else {
      simulation::material_xs.total      = 0.0;
      simulation::material_xs.absorption = 0.0;
      simulation::material_xs.fission    = 0.0;
      simulation::material_xs.nu_fission = 0.0;
    }

    // Find the distance to the nearest boundary
    double d_boundary;
    int surface_crossed;
    int lattice_translation[3];
    int next_level;
    distance_to_boundary(this, &d_boundary, &surface_crossed,
      lattice_translation, &next_level);

    // Sample a distance to collision
    double d_collision;
    if (type == static_cast<int>(ParticleType::electron) ||
        type == static_cast<int>(ParticleType::positron)) {
      d_collision = 0.0;
    } else if (simulation::material_xs.total == 0.0) {
      d_collision = INFINITY;
    } else {
      d_collision = -std::log(prn()) / simulation::material_xs.total;
    }

    // Select smaller of the two distances
    double distance = std::min(d_boundary, d_collision);

    // Advance particle
    for (int j = 0; j < n_coord; ++j) {
      // TODO: use Position
      coord[j].xyz[0] += distance * coord[j].uvw[0];
      coord[j].xyz[1] += distance * coord[j].uvw[1];
      coord[j].xyz[2] += distance * coord[j].uvw[2];
    }

    // Score track-length tallies
    if (!model::active_tracklength_tallies.empty()) {
      score_tracklength_tally(this, distance);
    }

    // Score track-length estimate of k-eff
    if (settings::run_mode == RUN_MODE_EIGENVALUE &&
        type == static_cast<int>(ParticleType::neutron)) {
      global_tally_tracklength += wgt * distance * simulation::material_xs.nu_fission;
    }

    // Score flux derivative accumulators for differential tallies.
    if (!model::active_tallies.empty()) {
      score_track_derivative(this, distance);
    }

    if (d_collision > d_boundary) {
      // ====================================================================
      // PARTICLE CROSSES SURFACE

      if (next_level > 0) n_coord = next_level;

      // Saving previous cell data
      for (int j = 0; j < n_coord; ++j) {
        last_cell[j] = coord[j].cell;
      }
      last_n_coord = n_coord;

      if (lattice_translation[0] != 0 || lattice_translation[1] != 0 ||
          lattice_translation[2] != 0) {
        // Particle crosses lattice boundary
        surface = ERROR_INT;
        cross_lattice(this, lattice_translation);
        event = EVENT_LATTICE;
      } else {
        // Particle crosses surface
        surface = surface_crossed;
        this->cross_surface();
        event = EVENT_SURFACE;
      }
      // Score cell to cell partial currents
      if (!model::active_surface_tallies.empty()) {
        score_surface_tally(this, model::active_surface_tallies);
      }
    } else {
      // ====================================================================
      // PARTICLE HAS COLLISION

      // Score collision estimate of keff
      if (settings::run_mode == RUN_MODE_EIGENVALUE &&
          type == static_cast<int>(ParticleType::neutron)) {
        global_tally_collision += wgt * simulation::material_xs.nu_fission
          / simulation::material_xs.total;
      }

      // Score surface current tallies -- this has to be done before the collision
      // since the direction of the particle will change and we need to use the
      // pre-collision direction to figure out what mesh surfaces were crossed

      if (!model::active_meshsurf_tallies.empty())
        score_surface_tally(this, model::active_meshsurf_tallies);

      // Clear surface component
      surface = ERROR_INT;

      if (settings::run_CE) {
        collision(this);
      } else {
        collision_mg(this);
      }

      // Score collision estimator tallies -- this is done after a collision
      // has occurred rather than before because we need information on the
      // outgoing energy for any tallies with an outgoing energy filter
      if (!model::active_collision_tallies.empty()) score_collision_tally(this);
      if (!model::active_analog_tallies.empty()) {
        if (settings::run_CE) {
          score_analog_tally_ce(this);
        } else {
          score_analog_tally_mg(this);
        }
      }

      // Reset banked weight during collision
      n_bank = 0;
      wgt_bank = 0.0;
      for (int& v : n_delayed_bank) v = 0;

      // Reset fission logical
      fission = false;

      // Save coordinates for tallying purposes
      std::copy(coord[0].xyz, coord[0].xyz + 3, last_xyz_current);

      // Set last material to none since cross sections will need to be
      // re-evaluated
      last_material = F90_NONE;

      // Set all uvws to base level -- right now, after a collision, only the
      // base level uvws are changed
      for (int j = 0; j < n_coord - 1; ++j) {
        if (coord[j + 1].rotated) {
          // If next level is rotated, apply rotation matrix
          const auto& m {model::cells[coord[j].cell]->rotation_};
          Direction u {coord[j].uvw};
          coord[j + 1].uvw[0] = m[3]*u.x + m[4]*u.y + m[5]*u.z;
          coord[j + 1].uvw[1] = m[6]*u.x + m[7]*u.y + m[8]*u.z;
          coord[j + 1].uvw[2] = m[9]*u.x + m[10]*u.y + m[11]*u.z;
        } else {
          // Otherwise, copy this level's direction
          std::copy(coord[j].uvw, coord[j].uvw + 3, coord[j + 1].uvw);
        }
      }

      // Score flux derivative accumulators for differential tallies.
      if (!model::active_tallies.empty()) score_collision_derivative(this);
    }

    // If particle has too many events, display warning and kill it
    ++n_event;
    if (n_event == MAX_EVENTS) {
      warning("Particle " + std::to_string(id) +
        " underwent maximum number of events.");
      alive = false;
    }

    // Check for secondary particles if this particle is dead
    if (!alive) {
      // If no secondary particles, break out of event loop
      if (n_secondary == 0) break;

      this->from_source(&secondary_bank[n_secondary - 1]);
      --n_secondary;
      n_event = 0;

      // Enter new particle in particle track file
      if (write_track) add_particle_track();
    }
  }

  // Finish particle track output.
  if (write_track) {
    write_particle_track(*this);
    finalize_particle_track(*this);
  }
}

void
Particle::cross_surface()
{
  int i_surface = std::abs(surface);
  // TODO: off-by-one
  const auto& surf {model::surfaces[i_surface - 1]};
  if (settings::verbosity >= 10 || simulation::trace) {
    write_message("    Crossing surface " + std::to_string(surf->id_));
  }

  if (surf->bc_ == BC_VACUUM && (settings::run_mode != RUN_MODE_PLOTTING)) {
    // =======================================================================
    // PARTICLE LEAKS OUT OF PROBLEM

    // Kill particle
    alive = false;

    // Score any surface current tallies -- note that the particle is moved
    // forward slightly so that if the mesh boundary is on the surface, it is
    // still processed

    if (!model::active_meshsurf_tallies.empty()) {
      // TODO: Find a better solution to score surface currents than
      // physically moving the particle forward slightly

      // TODO: Use Position
      coord[0].xyz[0] += TINY_BIT * coord[0].uvw[0];
      coord[0].xyz[1] += TINY_BIT * coord[0].uvw[1];
      coord[0].xyz[2] += TINY_BIT * coord[0].uvw[2];
      score_surface_tally(this, model::active_meshsurf_tallies);
    }

    // Score to global leakage tally
    global_tally_leakage += wgt;

    // Display message
    if (settings::verbosity >= 10 || simulation::trace) {
      write_message("    Leaked out of surface " + std::to_string(surf->id_));
    }
    return;

  } else if (surf->bc_ == BC_REFLECT && (settings::run_mode != RUN_MODE_PLOTTING)) {
    // =======================================================================
    // PARTICLE REFLECTS FROM SURFACE

    // Do not handle reflective boundary conditions on lower universes
    if (n_coord != 1) {
      this->mark_as_lost("Cannot reflect particle " + std::to_string(id) +
        " off surface in a lower universe.");
      return;
    }

    // Score surface currents since reflection causes the direction of the
    // particle to change. For surface filters, we need to score the tallies
    // twice, once before the particle's surface attribute has changed and
    // once after. For mesh surface filters, we need to artificially move
    // the particle slightly back in case the surface crossing is coincident
    // with a mesh boundary

    if (!model::active_surface_tallies.empty()) {
      score_surface_tally(this, model::active_surface_tallies);
    }


    if (!model::active_meshsurf_tallies.empty()) {
      Position r {coord[0].xyz};
      coord[0].xyz[0] -= TINY_BIT * coord[0].uvw[0];
      coord[0].xyz[1] -= TINY_BIT * coord[0].uvw[1];
      coord[0].xyz[2] -= TINY_BIT * coord[0].uvw[2];
      score_surface_tally(this, model::active_meshsurf_tallies);
      std::copy(&r.x, &r.x + 3, coord[0].xyz);
    }

    // Reflect particle off surface
    Direction u = surf->reflect(coord[0].xyz, coord[0].uvw);

    // Make sure new particle direction is normalized
    double norm = u.norm();
    coord[0].uvw[0] = u.x/norm;
    coord[0].uvw[1] = u.y/norm;
    coord[0].uvw[2] = u.z/norm;

    // Reassign particle's cell and surface
    coord[0].cell = last_cell[last_n_coord - 1];
    surface = -surface;

    // If a reflective surface is coincident with a lattice or universe
    // boundary, it is necessary to redetermine the particle's coordinates in
    // the lower universes.

    n_coord = 1;
    if (!find_cell(this, true)) {
      this->mark_as_lost("Couldn't find particle after reflecting from surface "
        + std::to_string(surf->id_) + ".");
      return;
    }

    // Set previous coordinate going slightly past surface crossing
    last_xyz_current[0] = coord[0].xyz[0] + TINY_BIT*coord[0].uvw[0];
    last_xyz_current[1] = coord[0].xyz[1] + TINY_BIT*coord[0].uvw[1];
    last_xyz_current[2] = coord[0].xyz[2] + TINY_BIT*coord[0].uvw[2];

    // Diagnostic message
    if (settings::verbosity >= 10 || simulation::trace) {
      write_message("    Reflected from surface " + std::to_string(surf->id_));
    }
    return;

  } else if (surf->bc_ == BC_PERIODIC && settings::run_mode != RUN_MODE_PLOTTING) {
    // =======================================================================
    // PERIODIC BOUNDARY

    // Do not handle periodic boundary conditions on lower universes
    if (n_coord != 1) {
      this->mark_as_lost("Cannot transfer particle " + std::to_string(id) +
        " across surface in a lower universe. Boundary conditions must be "
        "applied to root universe.");
      return;
    }

    // Score surface currents since reflection causes the direction of the
    // particle to change -- artificially move the particle slightly back in
    // case the surface crossing is coincident with a mesh boundary
    if (!model::active_meshsurf_tallies.empty()) {
      Position r {coord[0].xyz};
      coord[0].xyz[0] -= TINY_BIT * coord[0].uvw[0];
      coord[0].xyz[1] -= TINY_BIT * coord[0].uvw[1];
      coord[0].xyz[2] -= TINY_BIT * coord[0].uvw[2];
      score_surface_tally(this, model::active_meshsurf_tallies);
      std::copy(&r.x, &r.x + 3, coord[0].xyz);
    }

    // Get a pointer to the partner periodic surface
    auto surf_p = dynamic_cast<PeriodicSurface*>(surf);
    auto other = dynamic_cast<PeriodicSurface*>(
      model::surfaces[surf_p->i_periodic_]);

    // Adjust the particle's location and direction.
    Position r {coord[0].xyz};
    Direction u {coord[0].uvw};
    bool rotational = other->periodic_translate(surf_p, r, u);
    std::copy(&r.x, &r.x + 3, coord[0].xyz);
    std::copy(&u.x, &u.x + 3, coord[0].uvw);

    // Reassign particle's surface
    // TODO: off-by-one
    surface = rotational ?
      surf_p->i_periodic_ + 1 :
      std::copysign(surf_p->i_periodic_ + 1, surface);

    // Figure out what cell particle is in now
    n_coord = 1;

    if (!find_cell(this, true)) {
      this->mark_as_lost("Couldn't find particle after hitting periodic "
        "boundary on surface " + std::to_string(surf->id_) + ".");
      return;
    }

    // Set previous coordinate going slightly past surface crossing
    last_xyz_current[0] = coord[0].xyz[0] + TINY_BIT * coord[0].uvw[0];
    last_xyz_current[1] = coord[0].xyz[1] + TINY_BIT * coord[0].uvw[1];
    last_xyz_current[2] = coord[0].xyz[2] + TINY_BIT * coord[0].uvw[2];

    // Diagnostic message
    if (settings::verbosity >= 10 || simulation::trace) {
      write_message("    Hit periodic boundary on surface " +
        std::to_string(surf->id_));
    }
    return;
  }

  // ==========================================================================
  // SEARCH NEIGHBOR LISTS FOR NEXT CELL

#ifdef DAGMC
  if (settings::dagmc) {
    auto cellp = dynamic_cast<DAGCell*>(model::cells[last_cell[0]]);
    // TODO: off-by-one
    auto surfp = dynamic_cast<DAGSurface*>(model::surfaces[std::abs(surface) - 1]);
    int32_t i_cell = next_cell(cellp, surfp) - 1;
    // save material and temp
    last_material = material;
    last_sqrtkT = sqrtkT;
    // set new cell value
    coord[0].cell = i_cell;
    cell_instance = 0;
    // TODO: off-by-one
    int mat = model::cells[i_cell]->material_[0];
    material = (mat == MATERIAL_VOID) ? mat : mat + 1;
    sqrtkT = model::cells[i_cell]->sqrtkT_[0];
    return;
  }
#endif

  if (find_cell(this, true)) return;

  // ==========================================================================
  // COULDN'T FIND PARTICLE IN NEIGHBORING CELLS, SEARCH ALL CELLS

  // Remove lower coordinate levels and assignment of surface
  surface = ERROR_INT;
  n_coord = 1;
  bool found = find_cell(this, false);

  if (settings::run_mode != RUN_MODE_PLOTTING && (!found)) {
    // If a cell is still not found, there are two possible causes: 1) there is
    // a void in the model, and 2) the particle hit a surface at a tangent. If
    // the particle is really traveling tangent to a surface, if we move it
    // forward a tiny bit it should fix the problem.

    n_coord = 1;
    coord[0].xyz[0] += TINY_BIT * coord[0].uvw[0];
    coord[0].xyz[1] += TINY_BIT * coord[0].uvw[1];
    coord[0].xyz[2] += TINY_BIT * coord[0].uvw[2];

    // Couldn't find next cell anywhere! This probably means there is an actual
    // undefined region in the geometry.

    if (!find_cell(this, false)) {
      this->mark_as_lost("After particle " + std::to_string(id) +
        " crossed surface " + std::to_string(surf->id_) +
        " it could not be located in any cell and it did not leak.");
      return;
    }
  }
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
  simulation::n_lost_particles += 1;

  // Count the total number of simulated particles (on this processor)
  auto n = simulation::current_batch * settings::gen_per_batch * simulation::work;

  // Abort the simulation if the maximum number of lost particles has been
  // reached
  if (simulation::n_lost_particles >= MAX_LOST_PARTICLES &&
      simulation::n_lost_particles >= REL_MAX_LOST_PARTICLES*n) {
    fatal_error("Maximum number of lost particles has been reached.");
  }
}

void
Particle::write_restart() const
{
  // Dont write another restart file if in particle restart mode
  if (settings::run_mode == RUN_MODE_PARTICLE) return;

  // Set up file name
  std::stringstream filename;
  filename << settings::path_output << "particle_" << simulation::current_batch
    << '_' << id << ".h5";

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
    write_dataset(file_id, "current_batch", simulation::current_batch);
    write_dataset(file_id, "generations_per_batch", settings::gen_per_batch);
    write_dataset(file_id, "current_generation", simulation::current_gen);
    write_dataset(file_id, "n_particles", settings::n_particles);
    switch (settings::run_mode) {
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

    int64_t i = simulation::current_work;
    write_dataset(file_id, "weight", simulation::source_bank[i-1].wgt);
    write_dataset(file_id, "energy", simulation::source_bank[i-1].E);
    hsize_t dims[] {3};
    write_double(file_id, 1, dims, "xyz", simulation::source_bank[i-1].xyz, false);
    write_double(file_id, 1, dims, "uvw", simulation::source_bank[i-1].uvw, false);

    // Close file
    file_close(file_id);
  } // #pragma omp critical
}

//==============================================================================
// Fortran compatibility functions
//==============================================================================

void reset_coord(LocalCoord* c) { c->reset(); }
void particle_clear(Particle* p) { p->clear(); }
void particle_initialize(Particle* p) { p->initialize(); }
void particle_from_source(Particle* p, const Bank* src)
{
  p->from_source(src);
}
void particle_mark_as_lost(Particle* p, const char* message)
{
  p->mark_as_lost(message);
}
void particle_write_restart(Particle* p) { p->write_restart(); }

} // namespace openmc
