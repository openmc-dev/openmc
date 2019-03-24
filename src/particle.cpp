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
#include "openmc/photon.h"
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
  lattice = C_NONE;
  lattice_x = 0;
  lattice_y = 0;
  rotated = false;
}

//==============================================================================
// Particle implementation
//==============================================================================

Particle::Particle()
{
  // Create and clear coordinate levels
  coord_.resize(model::n_coord_levels);
  cell_last_.resize(model::n_coord_levels);
  clear();

  for (int& n : n_delayed_bank_) {
    n = 0;
  }

  // Create microscopic cross section caches
  neutron_xs_.resize(data::nuclides.size());
  photon_xs_.resize(data::elements.size());
}

void
Particle::clear()
{
  // reset any coordinate levels
  for (auto& level : coord_) level.reset();
  n_coord_ = 1;
}

void
Particle::create_secondary(Direction u, double E, Type type) const
{
  simulation::secondary_bank.emplace_back();

  auto& bank {simulation::secondary_bank.back()};
  bank.particle = type;
  bank.wgt = wgt_;
  bank.r = this->r();
  bank.u = u;
  bank.E = settings::run_CE ? E : g_;
}

void
Particle::from_source(const Bank* src)
{
  // reset some attributes
  this->clear();
  alive_ = true;
  surface_ = 0;
  cell_born_ = C_NONE;
  material_ = C_NONE;
  n_collision_ = 0;
  fission_ = false;

  // copy attributes from source bank site
  type_ = src->particle;
  wgt_ = src->wgt;
  wgt_last_ = src->wgt;
  this->r() = src->r;
  this->u() = src->u;
  r_last_current_ = src->r;
  r_last_ = src->r;
  u_last_ = src->u;
  if (settings::run_CE) {
    E_ = src->E;
    g_ = 0;
  } else {
    g_ = static_cast<int>(src->E);
    g_last_ = static_cast<int>(src->E);
    E_ = data::energy_bin_avg[g_ - 1];
  }
  E_last_ = E_;
}

void
Particle::transport()
{
  // Display message if high verbosity or trace is on
  if (settings::verbosity >= 9 || simulation::trace) {
     write_message("Simulating Particle " + std::to_string(id_));
  }

  // Initialize number of events to zero
  int n_event = 0;

  // Add paricle's starting weight to count for normalizing tallies later
  #pragma omp atomic
  simulation::total_weight += wgt_;

  // Force calculation of cross-sections by setting last energy to zero
  if (settings::run_CE) {
    for (auto& micro : neutron_xs_) micro.last_E = 0.0;
  }

  // Prepare to write out particle track.
  if (write_track_) add_particle_track();

  // Every particle starts with no accumulated flux derivative.
  if (!model::active_tallies.empty()) zero_flux_derivs();

  while (true) {
    // Set the random number stream
    if (type_ == Particle::Type::neutron) {
      prn_set_stream(STREAM_TRACKING);
    } else {
      prn_set_stream(STREAM_PHOTON);
    }

    // Store pre-collision particle properties
    wgt_last_ = wgt_;
    E_last_ = E_;
    u_last_ = this->u();
    r_last_ = this->r();

    // If the cell hasn't been determined based on the particle's location,
    // initiate a search for the current cell. This generally happens at the
    // beginning of the history and again for any secondary particles
    if (coord_[n_coord_ - 1].cell == C_NONE) {
      if (!find_cell(this, false)) {
        this->mark_as_lost("Could not find the cell containing particle "
          + std::to_string(id_));
        return;
      }

      // set birth cell attribute
      if (cell_born_ == C_NONE) cell_born_ = coord_[n_coord_ - 1].cell;
    }

    // Write particle track.
    if (write_track_) write_particle_track(*this);

    if (settings::check_overlaps) check_cell_overlap(this);

    // Calculate microscopic and macroscopic cross sections
    if (material_ != MATERIAL_VOID) {
      if (settings::run_CE) {
        if (material_ != material_last_ || sqrtkT_ != sqrtkT_last_) {
          // If the material is the same as the last material and the
          // temperature hasn't changed, we don't need to lookup cross
          // sections again.
          model::materials[material_]->calculate_xs(*this);
        }
      } else {
        // Get the MG data
        calculate_xs_c(material_, g_, sqrtkT_, this->u_local(),
          macro_xs_.total, macro_xs_.absorption, macro_xs_.nu_fission);

        // Finally, update the particle group while we have already checked
        // for if multi-group
        g_last_ = g_;
      }
    } else {
      macro_xs_.total      = 0.0;
      macro_xs_.absorption = 0.0;
      macro_xs_.fission    = 0.0;
      macro_xs_.nu_fission = 0.0;
    }

    // Find the distance to the nearest boundary
    auto boundary = distance_to_boundary(this);

    // Sample a distance to collision
    double d_collision;
    if (type_ == Particle::Type::electron ||
        type_ == Particle::Type::positron) {
      d_collision = 0.0;
    } else if (macro_xs_.total == 0.0) {
      d_collision = INFINITY;
    } else {
      d_collision = -std::log(prn()) / macro_xs_.total;
    }

    // Select smaller of the two distances
    double distance = std::min(boundary.distance, d_collision);

    // Advance particle
    for (int j = 0; j < n_coord_; ++j) {
      coord_[j].r += distance * coord_[j].u;
    }

    // Score track-length tallies
    if (!model::active_tracklength_tallies.empty()) {
      score_tracklength_tally(this, distance);
    }

    // Score track-length estimate of k-eff
    if (settings::run_mode == RUN_MODE_EIGENVALUE &&
        type_ == Particle::Type::neutron) {
      global_tally_tracklength += wgt_ * distance * macro_xs_.nu_fission;
    }

    // Score flux derivative accumulators for differential tallies.
    if (!model::active_tallies.empty()) {
      score_track_derivative(this, distance);
    }

    if (d_collision > boundary.distance) {
      // ====================================================================
      // PARTICLE CROSSES SURFACE

      // Set surface that particle is on and adjust coordinate levels
      surface_ = boundary.surface_index;
      n_coord_ = boundary.coord_level;

      // Saving previous cell data
      for (int j = 0; j < n_coord_; ++j) {
        cell_last_[j] = coord_[j].cell;
      }
      n_coord_last_ = n_coord_;

      if (boundary.lattice_translation[0] != 0 ||
          boundary.lattice_translation[1] != 0 ||
          boundary.lattice_translation[2] != 0) {
        // Particle crosses lattice boundary
        cross_lattice(this, boundary);
        event_ = EVENT_LATTICE;
      } else {
        // Particle crosses surface
        this->cross_surface();
        event_ = EVENT_SURFACE;
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
          type_ == Particle::Type::neutron) {
        global_tally_collision += wgt_ * macro_xs_.nu_fission
          / macro_xs_.total;
      }

      // Score surface current tallies -- this has to be done before the collision
      // since the direction of the particle will change and we need to use the
      // pre-collision direction to figure out what mesh surfaces were crossed

      if (!model::active_meshsurf_tallies.empty())
        score_surface_tally(this, model::active_meshsurf_tallies);

      // Clear surface component
      surface_ = 0;

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
      n_bank_ = 0;
      wgt_bank_ = 0.0;
      for (int& v : n_delayed_bank_) v = 0;

      // Reset fission logical
      fission_ = false;

      // Save coordinates for tallying purposes
      r_last_current_ = this->r();

      // Set last material to none since cross sections will need to be
      // re-evaluated
      material_last_ = C_NONE;

      // Set all directions to base level -- right now, after a collision, only
      // the base level directions are changed
      for (int j = 0; j < n_coord_ - 1; ++j) {
        if (coord_[j + 1].rotated) {
          // If next level is rotated, apply rotation matrix
          const auto& m {model::cells[coord_[j].cell]->rotation_};
          const auto& u {coord_[j].u};
          coord_[j + 1].u.x = m[3]*u.x + m[4]*u.y + m[5]*u.z;
          coord_[j + 1].u.y = m[6]*u.x + m[7]*u.y + m[8]*u.z;
          coord_[j + 1].u.z = m[9]*u.x + m[10]*u.y + m[11]*u.z;
        } else {
          // Otherwise, copy this level's direction
          coord_[j+1].u = coord_[j].u;
        }
      }

      // Score flux derivative accumulators for differential tallies.
      if (!model::active_tallies.empty()) score_collision_derivative(this);
    }

    // If particle has too many events, display warning and kill it
    ++n_event;
    if (n_event == MAX_EVENTS) {
      warning("Particle " + std::to_string(id_) +
        " underwent maximum number of events.");
      alive_ = false;
    }

    // Check for secondary particles if this particle is dead
    if (!alive_) {
      // If no secondary particles, break out of event loop
      if (simulation::secondary_bank.empty()) break;

      this->from_source(&simulation::secondary_bank.back());
      simulation::secondary_bank.pop_back();
      n_event = 0;

      // Enter new particle in particle track file
      if (write_track_) add_particle_track();
    }
  }

  // Finish particle track output.
  if (write_track_) {
    write_particle_track(*this);
    finalize_particle_track(*this);
  }
}

void
Particle::cross_surface()
{
  int i_surface = std::abs(surface_);
  // TODO: off-by-one
  const auto& surf {model::surfaces[i_surface - 1].get()};
  if (settings::verbosity >= 10 || simulation::trace) {
    write_message("    Crossing surface " + std::to_string(surf->id_));
  }

  if (surf->bc_ == BC_VACUUM && (settings::run_mode != RUN_MODE_PLOTTING)) {
    // =======================================================================
    // PARTICLE LEAKS OUT OF PROBLEM

    // Kill particle
    alive_ = false;

    // Score any surface current tallies -- note that the particle is moved
    // forward slightly so that if the mesh boundary is on the surface, it is
    // still processed

    if (!model::active_meshsurf_tallies.empty()) {
      // TODO: Find a better solution to score surface currents than
      // physically moving the particle forward slightly

      this->r() += TINY_BIT * this->u();
      score_surface_tally(this, model::active_meshsurf_tallies);
    }

    // Score to global leakage tally
    global_tally_leakage += wgt_;

    // Display message
    if (settings::verbosity >= 10 || simulation::trace) {
      write_message("    Leaked out of surface " + std::to_string(surf->id_));
    }
    return;

  } else if (surf->bc_ == BC_REFLECT && (settings::run_mode != RUN_MODE_PLOTTING)) {
    // =======================================================================
    // PARTICLE REFLECTS FROM SURFACE

    // Do not handle reflective boundary conditions on lower universes
    if (n_coord_ != 1) {
      this->mark_as_lost("Cannot reflect particle " + std::to_string(id_) +
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
      Position r {this->r()};
      this->r() -= TINY_BIT * this->u();
      score_surface_tally(this, model::active_meshsurf_tallies);
      this->r() = r;
    }

    // Reflect particle off surface
    Direction u = surf->reflect(this->r(), this->u());

    // Make sure new particle direction is normalized
    this->u() = u / u.norm();

    // Reassign particle's cell and surface
    coord_[0].cell = cell_last_[n_coord_last_ - 1];
    surface_ = -surface_;

    // If a reflective surface is coincident with a lattice or universe
    // boundary, it is necessary to redetermine the particle's coordinates in
    // the lower universes.

    n_coord_ = 1;
    if (!find_cell(this, true)) {
      this->mark_as_lost("Couldn't find particle after reflecting from surface "
        + std::to_string(surf->id_) + ".");
      return;
    }

    // Set previous coordinate going slightly past surface crossing
    r_last_current_ = this->r() + TINY_BIT*this->u();

    // Diagnostic message
    if (settings::verbosity >= 10 || simulation::trace) {
      write_message("    Reflected from surface " + std::to_string(surf->id_));
    }
    return;

  } else if (surf->bc_ == BC_PERIODIC && settings::run_mode != RUN_MODE_PLOTTING) {
    // =======================================================================
    // PERIODIC BOUNDARY

    // Do not handle periodic boundary conditions on lower universes
    if (n_coord_ != 1) {
      this->mark_as_lost("Cannot transfer particle " + std::to_string(id_) +
        " across surface in a lower universe. Boundary conditions must be "
        "applied to root universe.");
      return;
    }

    // Score surface currents since reflection causes the direction of the
    // particle to change -- artificially move the particle slightly back in
    // case the surface crossing is coincident with a mesh boundary
    if (!model::active_meshsurf_tallies.empty()) {
      Position r {this->r()};
      this->r() -= TINY_BIT * this->u();
      score_surface_tally(this, model::active_meshsurf_tallies);
      this->r() = r;
    }

    // Get a pointer to the partner periodic surface
    auto surf_p = dynamic_cast<PeriodicSurface*>(surf);
    auto other = dynamic_cast<PeriodicSurface*>(
      model::surfaces[surf_p->i_periodic_].get());

    // Adjust the particle's location and direction.
    bool rotational = other->periodic_translate(surf_p, this->r(), this->u());

    // Reassign particle's surface
    // TODO: off-by-one
    surface_ = rotational ?
      surf_p->i_periodic_ + 1 :
      std::copysign(surf_p->i_periodic_ + 1, surface_);

    // Figure out what cell particle is in now
    n_coord_ = 1;

    if (!find_cell(this, true)) {
      this->mark_as_lost("Couldn't find particle after hitting periodic "
        "boundary on surface " + std::to_string(surf->id_) + ".");
      return;
    }

    // Set previous coordinate going slightly past surface crossing
    r_last_current_ = this->r() + TINY_BIT*this->u();

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
    auto cellp = dynamic_cast<DAGCell*>(model::cells[cell_last_[0]].get());
    // TODO: off-by-one
    auto surfp = dynamic_cast<DAGSurface*>(model::surfaces[std::abs(surface_) - 1].get());
    int32_t i_cell = next_cell(cellp, surfp) - 1;
    // save material and temp
    material_last_ = material_;
    sqrtkT_last_ = sqrtkT_;
    // set new cell value
    coord_[0].cell = i_cell;
    cell_instance_ = 0;
    material_ = model::cells[i_cell]->material_[0];
    sqrtkT_ = model::cells[i_cell]->sqrtkT_[0];
    return;
  }
#endif

  if (find_cell(this, true)) return;

  // ==========================================================================
  // COULDN'T FIND PARTICLE IN NEIGHBORING CELLS, SEARCH ALL CELLS

  // Remove lower coordinate levels and assignment of surface
  surface_ = 0;
  n_coord_ = 1;
  bool found = find_cell(this, false);

  if (settings::run_mode != RUN_MODE_PLOTTING && (!found)) {
    // If a cell is still not found, there are two possible causes: 1) there is
    // a void in the model, and 2) the particle hit a surface at a tangent. If
    // the particle is really traveling tangent to a surface, if we move it
    // forward a tiny bit it should fix the problem.

    n_coord_ = 1;
    this->r() += TINY_BIT * this->u();

    // Couldn't find next cell anywhere! This probably means there is an actual
    // undefined region in the geometry.

    if (!find_cell(this, false)) {
      this->mark_as_lost("After particle " + std::to_string(id_) +
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
  alive_ = false;
  #pragma omp atomic
  simulation::n_lost_particles += 1;

  // Count the total number of simulated particles (on this processor)
  auto n = simulation::current_batch * settings::gen_per_batch *
    simulation::work_per_rank;

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
    << '_' << id_ << ".h5";

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
    write_dataset(file_id, "id", id_);
    write_dataset(file_id, "type", static_cast<int>(type_));

    int64_t i = simulation::current_work;
    write_dataset(file_id, "weight", simulation::source_bank[i-1].wgt);
    write_dataset(file_id, "energy", simulation::source_bank[i-1].E);
    write_dataset(file_id, "xyz", simulation::source_bank[i-1].r);
    write_dataset(file_id, "uvw", simulation::source_bank[i-1].u);

    // Close file
    file_close(file_id);
  } // #pragma omp critical
}

} // namespace openmc
