#include "openmc/particle.h"

#include <algorithm> // copy, min
#include <cmath>     // log, abs

#include <fmt/core.h>

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
#include "openmc/source.h"
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
LocalCoord::rotate(const std::vector<double>& rotation)
{
  this->r = this->r.rotate(rotation);
  this->u = this->u.rotate(rotation);
  this->rotated = true;
}

void
LocalCoord::reset()
{
  cell = C_NONE;
  universe = C_NONE;
  lattice = C_NONE;
  lattice_x = 0;
  lattice_y = 0;
  lattice_z = 0;
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
  // Reset any coordinate levels
  for (auto& level : coord_) level.reset();
  n_coord_ = 1;
}

void
Particle::create_secondary(double wgt, Direction u, double E, Type type)
{
  secondary_bank_.emplace_back();

  auto& bank {secondary_bank_.back()};
  bank.particle = type;
  bank.wgt = wgt;
  bank.r = this->r();
  bank.u = u;
  bank.E = settings::run_CE ? E : g_;

  n_bank_second_ += 1;
}

void
Particle::from_source(const Bank* src)
{
  // Reset some attributes
  this->clear();
  alive_ = true;
  surface_ = 0;
  cell_born_ = C_NONE;
  material_ = C_NONE;
  n_collision_ = 0;
  fission_ = false;
  std::fill(flux_derivs_.begin(), flux_derivs_.end(), 0.0);

  // Copy attributes from source bank site
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
    E_ = data::mg.energy_bin_avg_[g_];
  }
  E_last_ = E_;
}

void
Particle::event_calculate_xs()
{
  // Set the random number stream
  if (type_ == Particle::Type::neutron) {
    stream_ = STREAM_TRACKING;
  } else {
    stream_ = STREAM_PHOTON;
  }

  // Store pre-collision particle properties
  wgt_last_ = wgt_;
  E_last_ = E_;
  u_last_ = this->u();
  r_last_ = this->r();

  // Reset event variables
  event_ = TallyEvent::KILL;
  event_nuclide_ = NUCLIDE_NONE;
  event_mt_ = REACTION_NONE;

  // If the cell hasn't been determined based on the particle's location,
  // initiate a search for the current cell. This generally happens at the
  // beginning of the history and again for any secondary particles
  if (coord_[n_coord_ - 1].cell == C_NONE) {
    if (!exhaustive_find_cell(*this)) {
      this->mark_as_lost("Could not find the cell containing particle "
        + std::to_string(id_));
      return;
    }

    // Set birth cell attribute
    if (cell_born_ == C_NONE) cell_born_ = coord_[n_coord_ - 1].cell;
  }

  // Write particle track.
  if (write_track_) write_particle_track(*this);

  if (settings::check_overlaps) check_cell_overlap(*this);

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
      // Get the MG data; unlike the CE case above, we have to re-calculate
      // cross sections for every collision since the cross sections may
      // be angle-dependent
      data::mg.macro_xs_[material_].calculate_xs(*this);

      // Update the particle's group while we know we are multi-group
      g_last_ = g_;
    }
  } else {
    macro_xs_.total      = 0.0;
    macro_xs_.absorption = 0.0;
    macro_xs_.fission    = 0.0;
    macro_xs_.nu_fission = 0.0;
  }
}

void
Particle::event_advance()
{
  // Find the distance to the nearest boundary
  boundary_ = distance_to_boundary(*this);

  // Sample a distance to collision
  if (type_ == Particle::Type::electron ||
      type_ == Particle::Type::positron) {
    collision_distance_ = 0.0;
  } else if (macro_xs_.total == 0.0) {
    collision_distance_ = INFINITY;
  } else {
    collision_distance_ = -std::log(prn(this->current_seed())) / macro_xs_.total;
  }

  // Select smaller of the two distances
  double distance = std::min(boundary_.distance, collision_distance_);

  // Advance particle
  for (int j = 0; j < n_coord_; ++j) {
    coord_[j].r += distance * coord_[j].u;
  }

  // Score track-length tallies
  if (!model::active_tracklength_tallies.empty()) {
    score_tracklength_tally(*this, distance);
  }

  // Score track-length estimate of k-eff
  if (settings::run_mode == RunMode::EIGENVALUE &&
      type_ == Particle::Type::neutron) {
    keff_tally_tracklength_ += wgt_ * distance * macro_xs_.nu_fission;
  }

  // Score flux derivative accumulators for differential tallies.
  if (!model::active_tallies.empty()) {
    score_track_derivative(*this, distance);
  }
}

void
Particle::event_cross_surface()
{
  // Set surface that particle is on and adjust coordinate levels
  surface_ = boundary_.surface_index;
  n_coord_ = boundary_.coord_level;

  // Saving previous cell data
  for (int j = 0; j < n_coord_; ++j) {
    cell_last_[j] = coord_[j].cell;
  }
  n_coord_last_ = n_coord_;

  if (boundary_.lattice_translation[0] != 0 ||
      boundary_.lattice_translation[1] != 0 ||
      boundary_.lattice_translation[2] != 0) {
    // Particle crosses lattice boundary
    cross_lattice(*this, boundary_);
    event_ = TallyEvent::LATTICE;
  } else {
    // Particle crosses surface
    this->cross_surface();
    event_ = TallyEvent::SURFACE;
  }
  // Score cell to cell partial currents
  if (!model::active_surface_tallies.empty()) {
    score_surface_tally(*this, model::active_surface_tallies);
  }
}

void
Particle::event_collide()
{
  // Score collision estimate of keff
  if (settings::run_mode == RunMode::EIGENVALUE &&
      type_ == Particle::Type::neutron) {
    keff_tally_collision_ += wgt_ * macro_xs_.nu_fission
      / macro_xs_.total;
  }

  // Score surface current tallies -- this has to be done before the collision
  // since the direction of the particle will change and we need to use the
  // pre-collision direction to figure out what mesh surfaces were crossed

  if (!model::active_meshsurf_tallies.empty())
    score_surface_tally(*this, model::active_meshsurf_tallies);

  // Clear surface component
  surface_ = 0;

  if (settings::run_CE) {
    collision(*this);
  } else {
    collision_mg(*this);
  }

  // Score collision estimator tallies -- this is done after a collision
  // has occurred rather than before because we need information on the
  // outgoing energy for any tallies with an outgoing energy filter
  if (!model::active_collision_tallies.empty()) score_collision_tally(*this);
  if (!model::active_analog_tallies.empty()) {
    if (settings::run_CE) {
      score_analog_tally_ce(*this);
    } else {
      score_analog_tally_mg(*this);
    }
  }

  // Reset banked weight during collision
  n_bank_ = 0;
  n_bank_second_ = 0;
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
      coord_[j + 1].u = u.rotate(m);
    } else {
      // Otherwise, copy this level's direction
      coord_[j+1].u = coord_[j].u;
    }
  }

  // Score flux derivative accumulators for differential tallies.
  if (!model::active_tallies.empty()) score_collision_derivative(*this);
}

void
Particle::event_revive_from_secondary()
{
  // If particle has too many events, display warning and kill it
  ++n_event_;
  if (n_event_ == MAX_EVENTS) {
    warning("Particle " + std::to_string(id_) +
      " underwent maximum number of events.");
    alive_ = false;
  }

  // Check for secondary particles if this particle is dead
  if (!alive_) {
    // If no secondary particles, break out of event loop
    if (secondary_bank_.empty()) return;

    this->from_source(&secondary_bank_.back());
    secondary_bank_.pop_back();
    n_event_ = 0;

    // Enter new particle in particle track file
    if (write_track_) add_particle_track(*this);
  }
}

void
Particle::event_death()
{
  #ifdef DAGMC
  if (settings::dagmc) history_.reset();
  #endif

  // Finish particle track output.
  if (write_track_) {
    write_particle_track(*this);
    finalize_particle_track(*this);
  }

  // Contribute tally reduction variables to global accumulator
  #pragma omp atomic
  global_tally_absorption += keff_tally_absorption_;
  #pragma omp atomic
  global_tally_collision += keff_tally_collision_;
  #pragma omp atomic
  global_tally_tracklength += keff_tally_tracklength_;
  #pragma omp atomic
  global_tally_leakage += keff_tally_leakage_;

  // Reset particle tallies once accumulated
  keff_tally_absorption_  = 0.0;
  keff_tally_collision_   = 0.0;
  keff_tally_tracklength_ = 0.0;
  keff_tally_leakage_     = 0.0;

  // Record the number of progeny created by this particle.
  // This data will be used to efficiently sort the fission bank.
  if (settings::run_mode == RunMode::EIGENVALUE) {
    int64_t offset = id_ - 1 - simulation::work_index[mpi::rank];
    simulation::progeny_per_particle[offset] = n_progeny_;
  }
}


void
Particle::cross_surface()
{
  int i_surface = std::abs(surface_);
  // TODO: off-by-one
  const auto& surf {model::surfaces[i_surface - 1].get()};
  if (settings::verbosity >= 10 || trace_) {
    write_message(1, "    Crossing surface {}", surf->id_);
  }

  if (surf->surf_source_ && simulation::current_batch == settings::n_batches) {
    Particle::Bank site;
    site.r = this->r();
    site.u = this->u();
    site.E = this->E_;
    site.wgt = this->wgt_;
    site.delayed_group = this->delayed_group_;
    site.surf_id = surf->id_;
    site.particle = this->type_;
    site.parent_id = this->id_;
    site.progeny_id = this->n_progeny_;
    int64_t idx = simulation::surf_source_bank.thread_safe_append(site);
  }

  // Handle any applicable boundary conditions.
  if (surf->bc_ && settings::run_mode != RunMode::PLOTTING) {
    surf->bc_->handle_particle(*this, *surf);
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

  if (neighbor_list_find_cell(*this))
    return;

  // ==========================================================================
  // COULDN'T FIND PARTICLE IN NEIGHBORING CELLS, SEARCH ALL CELLS

  // Remove lower coordinate levels and assignment of surface
  surface_ = 0;
  n_coord_ = 1;
  bool found = exhaustive_find_cell(*this);

  if (settings::run_mode != RunMode::PLOTTING && (!found)) {
    // If a cell is still not found, there are two possible causes: 1) there is
    // a void in the model, and 2) the particle hit a surface at a tangent. If
    // the particle is really traveling tangent to a surface, if we move it
    // forward a tiny bit it should fix the problem.

    n_coord_ = 1;
    this->r() += TINY_BIT * this->u();

    // Couldn't find next cell anywhere! This probably means there is an actual
    // undefined region in the geometry.

    if (!exhaustive_find_cell(*this)) {
      this->mark_as_lost("After particle " + std::to_string(id_) +
        " crossed surface " + std::to_string(surf->id_) +
        " it could not be located in any cell and it did not leak.");
      return;
    }
  }
}

void
Particle::cross_vacuum_bc(const Surface& surf)
{
  // Kill the particle
  alive_ = false;

  // Score any surface current tallies -- note that the particle is moved
  // forward slightly so that if the mesh boundary is on the surface, it is
  // still processed

  if (!model::active_meshsurf_tallies.empty()) {
    // TODO: Find a better solution to score surface currents than
    // physically moving the particle forward slightly

    this->r() += TINY_BIT * this->u();
    score_surface_tally(*this, model::active_meshsurf_tallies);
  }

  // Score to global leakage tally
  keff_tally_leakage_ += wgt_;

  // Display message
  if (settings::verbosity >= 10 || trace_) {
    write_message(1, "    Leaked out of surface {}", surf.id_);
  }
}

void
Particle::cross_reflective_bc(const Surface& surf, Direction new_u)
{
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
    score_surface_tally(*this, model::active_surface_tallies);
  }

  if (!model::active_meshsurf_tallies.empty()) {
    Position r {this->r()};
    this->r() -= TINY_BIT * this->u();
    score_surface_tally(*this, model::active_meshsurf_tallies);
    this->r() = r;
  }

  // Set the new particle direction
  this->u() = new_u;

  // Reassign particle's cell and surface
  coord_[0].cell = cell_last_[n_coord_last_ - 1];
  surface_ = -surface_;

  // If a reflective surface is coincident with a lattice or universe
  // boundary, it is necessary to redetermine the particle's coordinates in
  // the lower universes.
  // (unless we're using a dagmc model, which has exactly one universe)
  if (!settings::dagmc) {
    n_coord_ = 1;
    if (!neighbor_list_find_cell(*this)) {
      this->mark_as_lost("Couldn't find particle after reflecting from surface "
                         + std::to_string(surf.id_) + ".");
      return;
    }
  }

  // Set previous coordinate going slightly past surface crossing
  r_last_current_ = this->r() + TINY_BIT*this->u();

  // Diagnostic message
  if (settings::verbosity >= 10 || trace_) {
    write_message(1, "    Reflected from surface {}", surf.id_);
  }
}

void
Particle::cross_periodic_bc(const Surface& surf, Position new_r,
                            Direction new_u, int new_surface)
{
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
    score_surface_tally(*this, model::active_meshsurf_tallies);
    this->r() = r;
  }

  // Adjust the particle's location and direction.
  r() = new_r;
  u() = new_u;

  // Reassign particle's surface
  surface_ = new_surface;

  // Figure out what cell particle is in now
  n_coord_ = 1;

  if (!neighbor_list_find_cell(*this)) {
    this->mark_as_lost("Couldn't find particle after hitting periodic "
      "boundary on surface " + std::to_string(surf.id_) + ". The normal vector "
      "of one periodic surface may need to be reversed.");
    return;
  }

  // Set previous coordinate going slightly past surface crossing
  r_last_current_ = this->r() + TINY_BIT*this->u();

  // Diagnostic message
  if (settings::verbosity >= 10 || trace_) {
    write_message(1, "    Hit periodic boundary on surface {}", surf.id_);
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
  if (simulation::n_lost_particles >= settings::max_lost_particles &&
      simulation::n_lost_particles >= settings::rel_max_lost_particles*n) {
    fatal_error("Maximum number of lost particles has been reached.");
  }
}

void
Particle::write_restart() const
{
  // Dont write another restart file if in particle restart mode
  if (settings::run_mode == RunMode::PARTICLE) return;

  // Set up file name
  auto filename = fmt::format("{}particle_{}_{}.h5", settings::path_output,
    simulation::current_batch, id_);

  #pragma omp critical (WriteParticleRestart)
  {
    // Create file
    hid_t file_id = file_open(filename, 'w');

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
      case RunMode::FIXED_SOURCE:
        write_dataset(file_id, "run_mode", "fixed source");
        break;
      case RunMode::EIGENVALUE:
        write_dataset(file_id, "run_mode", "eigenvalue");
        break;
      case RunMode::PARTICLE:
        write_dataset(file_id, "run_mode", "particle restart");
        break;
      default:
        break;
    }
    write_dataset(file_id, "id", id_);
    write_dataset(file_id, "type", static_cast<int>(type_));

    int64_t i = current_work_;
    if (settings::run_mode == RunMode::EIGENVALUE) {
      //take source data from primary bank for eigenvalue simulation
      write_dataset(file_id, "weight", simulation::source_bank[i-1].wgt);
      write_dataset(file_id, "energy", simulation::source_bank[i-1].E);
      write_dataset(file_id, "xyz", simulation::source_bank[i-1].r);
      write_dataset(file_id, "uvw", simulation::source_bank[i-1].u);
    } else if (settings::run_mode == RunMode::FIXED_SOURCE) {
      // re-sample using rng random number seed used to generate source particle
      int64_t id = (simulation::total_gen + overall_generation() - 1)*settings::n_particles +
        simulation::work_index[mpi::rank] + i;
      uint64_t seed = init_seed(id, STREAM_SOURCE);
      // re-sample source site
      auto site = sample_external_source(&seed);
      write_dataset(file_id, "weight", site.wgt);
      write_dataset(file_id, "energy", site.E);
      write_dataset(file_id, "xyz", site.r);
      write_dataset(file_id, "uvw", site.u);
    }

    // Close file
    file_close(file_id);
  } // #pragma omp critical
}

std::string particle_type_to_str(Particle::Type type)
{
  switch (type) {
    case Particle::Type::neutron:
      return "neutron";
    case Particle::Type::photon:
      return "photon";
    case Particle::Type::electron:
      return "electron";
    case Particle::Type::positron:
      return "positron";
  }
  UNREACHABLE();
}

Particle::Type str_to_particle_type(std::string str)
{
  if (str == "neutron") {
    return Particle::Type::neutron;
  } else if (str == "photon") {
    return Particle::Type::photon;
  } else if (str == "electron") {
    return Particle::Type::electron;
  } else if (str == "positron") {
    return Particle::Type::positron;
  } else {
    throw std::invalid_argument{fmt::format("Invalid particle name: {}", str)};
  }
}

} // namespace openmc
