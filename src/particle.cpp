#include "openmc/particle.h"

#include <algorithm> // copy, min
#include <cmath>     // log, abs

#include <fmt/core.h>

#include "openmc/bank.h"
#include "openmc/capi.h"
#include "openmc/cell.h"
#include "openmc/constants.h"
#include "openmc/dagmc.h"
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
#include "openmc/simulation.h"
#include "openmc/source.h"
#include "openmc/surface.h"
#include "openmc/tallies/derivative.h"
#include "openmc/tallies/tally.h"
#include "openmc/tallies/tally_scoring.h"
#include "openmc/track_output.h"

#ifdef DAGMC
#include "DagMC.hpp"
#endif

namespace openmc {

double Particle::speed() const
{
  // Determine mass in eV/c^2
  double mass;
  switch (this->type()) {
  case ParticleType::neutron:
    mass = MASS_NEUTRON_EV;
    break;
  case ParticleType::photon:
    mass = 0.0;
    break;
  case ParticleType::electron:
  case ParticleType::positron:
    mass = MASS_ELECTRON_EV;
    break;
  }

  // Calculate inverse of Lorentz factor
  const double inv_gamma = mass / (this->E() + mass);

  // Calculate speed via v = c * sqrt(1 - γ^-2)
  return C_LIGHT * std::sqrt(1 - inv_gamma * inv_gamma);
}

void Particle::create_secondary(
  double wgt, Direction u, double E, ParticleType type)
{
  // If energy is below cutoff for this particle, don't create secondary
  // particle
  if (E < settings::energy_cutoff[static_cast<int>(type)]) {
    return;
  }

  secondary_bank().emplace_back();

  auto& bank {secondary_bank().back()};
  bank.particle = type;
  bank.wgt = wgt;
  bank.r = r();
  bank.u = u;
  bank.E = settings::run_CE ? E : g();
  bank.time = time();

  n_bank_second() += 1;
}

void Particle::from_source(const SourceSite* src)
{
  // Reset some attributes
  clear();
  alive() = true;
  surface() = 0;
  cell_born() = C_NONE;
  material() = C_NONE;
  n_collision() = 0;
  fission() = false;
  zero_flux_derivs();

  // Copy attributes from source bank site
  type() = src->particle;
  wgt() = src->wgt;
  wgt_last() = src->wgt;
  r() = src->r;
  u() = src->u;
  r_last_current() = src->r;
  r_last() = src->r;
  u_last() = src->u;
  if (settings::run_CE) {
    E() = src->E;
    g() = 0;
  } else {
    g() = static_cast<int>(src->E);
    g_last() = static_cast<int>(src->E);
    E() = data::mg.energy_bin_avg_[g()];
  }
  E_last() = E();
  time() = src->time;
  time_last() = src->time;
}

void Particle::event_calculate_xs()
{
  // Set the random number stream
  stream() = STREAM_TRACKING;

  // Store pre-collision particle properties
  wgt_last() = wgt();
  E_last() = E();
  u_last() = u();
  r_last() = r();
  time_last() = time();

  // Reset event variables
  event() = TallyEvent::KILL;
  event_nuclide() = NUCLIDE_NONE;
  event_mt() = REACTION_NONE;

  // If the cell hasn't been determined based on the particle's location,
  // initiate a search for the current cell. This generally happens at the
  // beginning of the history and again for any secondary particles
  if (coord(n_coord() - 1).cell == C_NONE) {
    if (!exhaustive_find_cell(*this)) {
      mark_as_lost(
        "Could not find the cell containing particle " + std::to_string(id()));
      return;
    }

    // Set birth cell attribute
    if (cell_born() == C_NONE)
      cell_born() = coord(n_coord() - 1).cell;
  }

  // Write particle track.
  if (write_track())
    write_particle_track(*this);

  if (settings::check_overlaps)
    check_cell_overlap(*this);

  // Calculate microscopic and macroscopic cross sections
  if (material() != MATERIAL_VOID) {
    if (settings::run_CE) {
      if (material() != material_last() || sqrtkT() != sqrtkT_last()) {
        // If the material is the same as the last material and the
        // temperature hasn't changed, we don't need to lookup cross
        // sections again.
        model::materials[material()]->calculate_xs(*this);
      }
    } else {
      // Get the MG data; unlike the CE case above, we have to re-calculate
      // cross sections for every collision since the cross sections may
      // be angle-dependent
      data::mg.macro_xs_[material()].calculate_xs(*this);

      // Update the particle's group while we know we are multi-group
      g_last() = g();
    }
  } else {
    macro_xs().total = 0.0;
    macro_xs().absorption = 0.0;
    macro_xs().fission = 0.0;
    macro_xs().nu_fission = 0.0;
  }
}

void Particle::event_advance()
{
  // Find the distance to the nearest boundary

  // Sample a distance to collision
  if (type() == ParticleType::electron || type() == ParticleType::positron) {
    collision_distance() = 0.0;
  } else if (macro_xs().total == 0.0) {
    collision_distance() = INFINITY;
  } else {
    collision_distance() = -std::log(prn(current_seed())) / macro_xs().total;
  }

  boundary() = distance_to_boundary(*this);

  // Select smaller of the two distances
  double distance = std::min(boundary().distance, collision_distance());

  // Advance particle in space and time
  for (int j = 0; j < n_coord(); ++j) {
    coord(j).r += distance * coord(j).u;
  }
  this->time() += distance / this->speed();

  // Score track-length tallies
  if (!model::active_tracklength_tallies.empty()) {
    score_tracklength_tally(*this, distance);
  }

  // Score track-length estimate of k-eff
  if (settings::run_mode == RunMode::EIGENVALUE &&
      type() == ParticleType::neutron) {
    keff_tally_tracklength() += wgt() * distance * macro_xs().nu_fission;
  }

  // Score flux derivative accumulators for differential tallies.
  if (!model::active_tallies.empty()) {
    score_track_derivative(*this, distance);
  }
}

void Particle::event_cross_surface()
{
  // Set surface that particle is on and adjust coordinate levels
  surface() = boundary().surface_index;
  n_coord() = boundary().coord_level;

  // Saving previous cell data
  for (int j = 0; j < n_coord(); ++j) {
    cell_last(j) = coord(j).cell;
  }
  n_coord_last() = n_coord();

  if (boundary().lattice_translation[0] != 0 ||
      boundary().lattice_translation[1] != 0 ||
      boundary().lattice_translation[2] != 0) {
    // Particle crosses lattice boundary
    cross_lattice(*this, boundary());
    event() = TallyEvent::LATTICE;
  } else {
    // Particle crosses surface
    cross_surface();
    event() = TallyEvent::SURFACE;
  }
  // Score cell to cell partial currents
  if (!model::active_surface_tallies.empty()) {
    score_surface_tally(*this, model::active_surface_tallies);
  }
}

void Particle::event_collide()
{
  // Score collision estimate of keff
  if (settings::run_mode == RunMode::EIGENVALUE &&
      type() == ParticleType::neutron) {
    keff_tally_collision() += wgt() * macro_xs().nu_fission / macro_xs().total;
  }

  // Score surface current tallies -- this has to be done before the collision
  // since the direction of the particle will change and we need to use the
  // pre-collision direction to figure out what mesh surfaces were crossed

  if (!model::active_meshsurf_tallies.empty())
    score_surface_tally(*this, model::active_meshsurf_tallies);

  // Clear surface component
  surface() = 0;

  if (settings::run_CE) {
    collision(*this);
  } else {
    collision_mg(*this);
  }

  // Score collision estimator tallies -- this is done after a collision
  // has occurred rather than before because we need information on the
  // outgoing energy for any tallies with an outgoing energy filter
  if (!model::active_collision_tallies.empty())
    score_collision_tally(*this);
  if (!model::active_analog_tallies.empty()) {
    if (settings::run_CE) {
      score_analog_tally_ce(*this);
    } else {
      score_analog_tally_mg(*this);
    }
  }

  // Reset banked weight during collision
  n_bank() = 0;
  n_bank_second() = 0;
  wgt_bank() = 0.0;
  zero_delayed_bank();

  // Reset fission logical
  fission() = false;

  // Save coordinates for tallying purposes
  r_last_current() = r();

  // Set last material to none since cross sections will need to be
  // re-evaluated
  material_last() = C_NONE;

  // Set all directions to base level -- right now, after a collision, only
  // the base level directions are changed
  for (int j = 0; j < n_coord() - 1; ++j) {
    if (coord(j + 1).rotated) {
      // If next level is rotated, apply rotation matrix
      const auto& m {model::cells[coord(j).cell]->rotation_};
      const auto& u {coord(j).u};
      coord(j + 1).u = u.rotate(m);
    } else {
      // Otherwise, copy this level's direction
      coord(j + 1).u = coord(j).u;
    }
  }

  // Score flux derivative accumulators for differential tallies.
  if (!model::active_tallies.empty())
    score_collision_derivative(*this);

#ifdef DAGMC
  history().reset();
#endif
}

void Particle::event_revive_from_secondary()
{
  // If particle has too many events, display warning and kill it
  ++n_event();
  if (n_event() == MAX_EVENTS) {
    warning("Particle " + std::to_string(id()) +
            " underwent maximum number of events.");
    alive() = false;
  }

  // Check for secondary particles if this particle is dead
  if (!alive()) {
    // If no secondary particles, break out of event loop
    if (secondary_bank().empty())
      return;

    from_source(&secondary_bank().back());
    secondary_bank().pop_back();
    n_event() = 0;

    // Enter new particle in particle track file
    if (write_track())
      add_particle_track(*this);
  }
}

void Particle::event_death()
{
#ifdef DAGMC
  history().reset();
#endif

  // Finish particle track output.
  if (write_track()) {
    write_particle_track(*this);
    finalize_particle_track(*this);
  }

// Contribute tally reduction variables to global accumulator
#pragma omp atomic
  global_tally_absorption += keff_tally_absorption();
#pragma omp atomic
  global_tally_collision += keff_tally_collision();
#pragma omp atomic
  global_tally_tracklength += keff_tally_tracklength();
#pragma omp atomic
  global_tally_leakage += keff_tally_leakage();

  // Reset particle tallies once accumulated
  keff_tally_absorption() = 0.0;
  keff_tally_collision() = 0.0;
  keff_tally_tracklength() = 0.0;
  keff_tally_leakage() = 0.0;

  // Record the number of progeny created by this particle.
  // This data will be used to efficiently sort the fission bank.
  if (settings::run_mode == RunMode::EIGENVALUE) {
    int64_t offset = id() - 1 - simulation::work_index[mpi::rank];
    simulation::progeny_per_particle[offset] = n_progeny();
  }
}

void Particle::cross_surface()
{
  int i_surface = std::abs(surface());
  // TODO: off-by-one
  const auto& surf {model::surfaces[i_surface - 1].get()};
  if (settings::verbosity >= 10 || trace()) {
    write_message(1, "    Crossing surface {}", surf->id_);
  }

  if (surf->surf_source_ && simulation::current_batch == settings::n_batches) {
    SourceSite site;
    site.r = r();
    site.u = u();
    site.E = E();
    site.time = time();
    site.wgt = wgt();
    site.delayed_group = delayed_group();
    site.surf_id = surf->id_;
    site.particle = type();
    site.parent_id = id();
    site.progeny_id = n_progeny();
    int64_t idx = simulation::surf_source_bank.thread_safe_append(site);
  }

// if we're crossing a CSG surface, make sure the DAG history is reset
#ifdef DAGMC
  if (surf->geom_type_ == GeometryType::CSG)
    history().reset();
#endif

  // Handle any applicable boundary conditions.
  if (surf->bc_ && settings::run_mode != RunMode::PLOTTING) {
    surf->bc_->handle_particle(*this, *surf);
    return;
  }

  // ==========================================================================
  // SEARCH NEIGHBOR LISTS FOR NEXT CELL

#ifdef DAGMC
  // in DAGMC, we know what the next cell should be
  if (surf->geom_type_ == GeometryType::DAG) {
    auto surfp = dynamic_cast<DAGSurface*>(surf);
    auto cellp =
      dynamic_cast<DAGCell*>(model::cells[cell_last(n_coord() - 1)].get());
    auto univp = static_cast<DAGUniverse*>(
      model::universes[coord(n_coord() - 1).universe].get());
    // determine the next cell for this crossing
    int32_t i_cell = next_cell(univp, cellp, surfp) - 1;
    // save material and temp
    material_last() = material();
    sqrtkT_last() = sqrtkT();
    // set new cell value
    coord(n_coord() - 1).cell = i_cell;
    cell_instance() = 0;
    material() = model::cells[i_cell]->material_[0];
    sqrtkT() = model::cells[i_cell]->sqrtkT_[0];
    return;
  }
#endif

  if (surf->is_triso_surface_) {
    if (surface() > 0){
      for (int i = n_coord(); i < model::n_coord_levels; i++) {
        coord(i).reset();
      }
      coord(n_coord() - 1).cell = model::cell_map[model::surfaces[i_surface - 1]->triso_base_index_];
    } else if (surface() < 0) {
      for (int i = n_coord(); i < model::n_coord_levels; i++) {
        coord(i).reset();
      }
      if (model::surfaces[i_surface - 1]->triso_particle_index_==-1) {
        fatal_error(fmt::format("Particle cell of surface {} is not defined", model::surfaces[i_surface - 1]->id_));
      }
      coord(n_coord() - 1).cell = model::cell_map[model::surfaces[i_surface - 1]->triso_particle_index_];
    }

    //find material
    bool found=true;
    int i_cell = coord(n_coord() - 1).cell;
    for (;; ++n_coord()) {
      if (i_cell == C_NONE) {
        int i_universe = coord(n_coord() - 1).universe;
        const auto& univ {model::universes[i_universe]};

        if (univ->filled_with_triso_base_ != -1) {
          coord(n_coord() - 1).cell = model::cell_map[univ->filled_with_triso_base_];
          found=true;
        } else {
          found = univ->find_cell(*this);
        }
        if (!found) {
          break;
        }
      }

      i_cell = coord(n_coord() - 1).cell;

      Cell& c {*model::cells[i_cell]};
      if (c.type_ == Fill::MATERIAL) {
        // Found a material cell which means this is the lowest coord level.

        cell_instance() = 0;
        // Find the distribcell instance number.
        if (c.distribcell_index_ >= 0) {
          cell_instance() = cell_instance_at_level(*this, n_coord() - 1);
        }

        // Set the material and temperature.
        material_last() = material();
        if (c.material_.size() > 1) {
          material() = c.material_[cell_instance()];
        } else {
          material() = c.material_[0];
        }
        sqrtkT_last() = sqrtkT();
        if (c.sqrtkT_.size() > 1) {
          sqrtkT() = c.sqrtkT_[cell_instance()];
        } else {
          sqrtkT() = c.sqrtkT_[0];
        }
        return;

      } else if (c.type_ == Fill::UNIVERSE) {
        //========================================================================
        //! Found a lower universe, update this coord level then search the next.

        // Set the lower coordinate level universe.
        auto& coor {coord(n_coord())};
        coor.universe = c.fill_;

        // Set the position and direction.
        coor.r = r_local();
        coor.u = u_local();

        // Apply translation.
        coor.r -= c.translation_;

        // Apply rotation.
        if (!c.rotation_.empty()) {
          coor.rotate(c.rotation_);
        }
        i_cell = C_NONE;
      }
    }
  } else {
      if (neighbor_list_find_cell(*this))
      return;
  }

  // ==========================================================================
  // COULDN'T FIND PARTICLE IN NEIGHBORING CELLS, SEARCH ALL CELLS

  // Remove lower coordinate levels and assignment of surface
  surface() = 0;
  n_coord() = 1;
  bool found = exhaustive_find_cell(*this);

  if (settings::run_mode != RunMode::PLOTTING && (!found)) {
    // If a cell is still not found, there are two possible causes: 1) there is
    // a void in the model, and 2) the particle hit a surface at a tangent. If
    // the particle is really traveling tangent to a surface, if we move it
    // forward a tiny bit it should fix the problem.

    n_coord() = 1;
    r() += TINY_BIT * u();

    // Couldn't find next cell anywhere! This probably means there is an actual
    // undefined region in the geometry.

    if (!exhaustive_find_cell(*this)) {
      mark_as_lost("After particle " + std::to_string(id()) +
                   " crossed surface " + std::to_string(surf->id_) +
                   " it could not be located in any cell and it did not leak.");
      return;
    }
  }
}

void Particle::cross_vacuum_bc(const Surface& surf)
{
  // Kill the particle
  alive() = false;

  // Score any surface current tallies -- note that the particle is moved
  // forward slightly so that if the mesh boundary is on the surface, it is
  // still processed

  if (!model::active_meshsurf_tallies.empty()) {
    // TODO: Find a better solution to score surface currents than
    // physically moving the particle forward slightly

    r() += TINY_BIT * u();
    score_surface_tally(*this, model::active_meshsurf_tallies);
  }

  // Score to global leakage tally
  keff_tally_leakage() += wgt();

  // Display message
  if (settings::verbosity >= 10 || trace()) {
    write_message(1, "    Leaked out of surface {}", surf.id_);
  }
}

void Particle::cross_reflective_bc(const Surface& surf, Direction new_u)
{
  // Do not handle reflective boundary conditions on lower universes
  if (n_coord() != 1) {
    mark_as_lost("Cannot reflect particle " + std::to_string(id()) +
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
    this->r() -= TINY_BIT * u();
    score_surface_tally(*this, model::active_meshsurf_tallies);
    this->r() = r;
  }

  // Set the new particle direction
  u() = new_u;

  // Reassign particle's cell and surface
  coord(0).cell = cell_last(n_coord_last() - 1);
  surface() = -surface();

  // If a reflective surface is coincident with a lattice or universe
  // boundary, it is necessary to redetermine the particle's coordinates in
  // the lower universes.
  // (unless we're using a dagmc model, which has exactly one universe)
  n_coord() = 1;
  if (surf.geom_type_ != GeometryType::DAG && !neighbor_list_find_cell(*this)) {
    this->mark_as_lost("Couldn't find particle after reflecting from surface " +
                       std::to_string(surf.id_) + ".");
    return;
  }

  // Set previous coordinate going slightly past surface crossing
  r_last_current() = r() + TINY_BIT * u();

  // Diagnostic message
  if (settings::verbosity >= 10 || trace()) {
    write_message(1, "    Reflected from surface {}", surf.id_);
  }
}

void Particle::cross_periodic_bc(
  const Surface& surf, Position new_r, Direction new_u, int new_surface)
{
  // Do not handle periodic boundary conditions on lower universes
  if (n_coord() != 1) {
    mark_as_lost(
      "Cannot transfer particle " + std::to_string(id()) +
      " across surface in a lower universe. Boundary conditions must be "
      "applied to root universe.");
    return;
  }

  // Score surface currents since reflection causes the direction of the
  // particle to change -- artificially move the particle slightly back in
  // case the surface crossing is coincident with a mesh boundary
  if (!model::active_meshsurf_tallies.empty()) {
    Position r {this->r()};
    this->r() -= TINY_BIT * u();
    score_surface_tally(*this, model::active_meshsurf_tallies);
    this->r() = r;
  }

  // Adjust the particle's location and direction.
  r() = new_r;
  u() = new_u;

  // Reassign particle's surface
  surface() = new_surface;

  // Figure out what cell particle is in now
  n_coord() = 1;

  if (!neighbor_list_find_cell(*this)) {
    mark_as_lost("Couldn't find particle after hitting periodic "
                 "boundary on surface " +
                 std::to_string(surf.id_) +
                 ". The normal vector "
                 "of one periodic surface may need to be reversed.");
    return;
  }

  // Set previous coordinate going slightly past surface crossing
  r_last_current() = r() + TINY_BIT * u();

  // Diagnostic message
  if (settings::verbosity >= 10 || trace()) {
    write_message(1, "    Hit periodic boundary on surface {}", surf.id_);
  }
}

void Particle::mark_as_lost(const char* message)
{
  // Print warning and write lost particle file
  warning(message);
  write_restart();

  // Increment number of lost particles
  alive() = false;
#pragma omp atomic
  simulation::n_lost_particles += 1;

  // Count the total number of simulated particles (on this processor)
  auto n = simulation::current_batch * settings::gen_per_batch *
           simulation::work_per_rank;

  // Abort the simulation if the maximum number of lost particles has been
  // reached
  if (simulation::n_lost_particles >= settings::max_lost_particles &&
      simulation::n_lost_particles >= settings::rel_max_lost_particles * n) {
    fatal_error("Maximum number of lost particles has been reached.");
  }
}

void Particle::write_restart() const
{
  // Dont write another restart file if in particle restart mode
  if (settings::run_mode == RunMode::PARTICLE)
    return;

  // Set up file name
  auto filename = fmt::format("{}particle_{}_{}.h5", settings::path_output,
    simulation::current_batch, id());

#pragma omp critical(WriteParticleRestart)
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
    write_dataset(file_id, "id", id());
    write_dataset(file_id, "type", static_cast<int>(type()));

    int64_t i = current_work();
    if (settings::run_mode == RunMode::EIGENVALUE) {
      // take source data from primary bank for eigenvalue simulation
      write_dataset(file_id, "weight", simulation::source_bank[i - 1].wgt);
      write_dataset(file_id, "energy", simulation::source_bank[i - 1].E);
      write_dataset(file_id, "xyz", simulation::source_bank[i - 1].r);
      write_dataset(file_id, "uvw", simulation::source_bank[i - 1].u);
      write_dataset(file_id, "time", simulation::source_bank[i - 1].time);
    } else if (settings::run_mode == RunMode::FIXED_SOURCE) {
      // re-sample using rng random number seed used to generate source particle
      int64_t id = (simulation::total_gen + overall_generation() - 1) *
                     settings::n_particles +
                   simulation::work_index[mpi::rank] + i;
      uint64_t seed = init_seed(id, STREAM_SOURCE);
      // re-sample source site
      auto site = sample_external_source(&seed);
      write_dataset(file_id, "weight", site.wgt);
      write_dataset(file_id, "energy", site.E);
      write_dataset(file_id, "xyz", site.r);
      write_dataset(file_id, "uvw", site.u);
      write_dataset(file_id, "time", site.time);
    }

    // Close file
    file_close(file_id);
  } // #pragma omp critical
}

std::string particle_type_to_str(ParticleType type)
{
  switch (type) {
  case ParticleType::neutron:
    return "neutron";
  case ParticleType::photon:
    return "photon";
  case ParticleType::electron:
    return "electron";
  case ParticleType::positron:
    return "positron";
  }
  UNREACHABLE();
}

ParticleType str_to_particle_type(std::string str)
{
  if (str == "neutron") {
    return ParticleType::neutron;
  } else if (str == "photon") {
    return ParticleType::photon;
  } else if (str == "electron") {
    return ParticleType::electron;
  } else if (str == "positron") {
    return ParticleType::positron;
  } else {
    throw std::invalid_argument {fmt::format("Invalid particle name: {}", str)};
  }
}

} // namespace openmc
