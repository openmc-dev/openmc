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
#include "openmc/particle_data.h"
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
#include "openmc/weight_windows.h"

#ifdef OPENMC_DAGMC_ENABLED
#include "DagMC.hpp"
#endif

namespace openmc {

//==============================================================================
// Particle implementation
//==============================================================================

double Particle::speed() const
{
  if (settings::run_CE) {
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
    // Equivalent to C * sqrt(1-(m/(m+E))^2) without problem at E<<m:
    return C_LIGHT * std::sqrt(this->E() * (this->E() + 2 * mass)) /
           (this->E() + mass);
  } else {
    auto& macro_xs = data::mg.macro_xs_[this->material()];
    int macro_t = this->mg_xs_cache().t;
    int macro_a = macro_xs.get_angle_index(this->u());
    return 1.0 / macro_xs.get_xs(MgxsType::INVERSE_VELOCITY, this->g(), nullptr,
                   nullptr, nullptr, macro_t, macro_a);
  }
}

bool Particle::create_secondary(
  double wgt, Direction u, double E, ParticleType type)
{
  // If energy is below cutoff for this particle, don't create secondary
  // particle
  if (E < settings::energy_cutoff[static_cast<int>(type)]) {
    return false;
  }

  auto& bank = secondary_bank().emplace_back();
  bank.particle = type;
  bank.wgt = wgt;
  bank.r = r();
  bank.u = u;
  bank.E = settings::run_CE ? E : g();
  bank.time = time();
  bank_second_E() += bank.E;
  return true;
}

void Particle::split(double wgt)
{
  auto& bank = secondary_bank().emplace_back();
  bank.particle = type();
  bank.wgt = wgt;
  bank.r = r();
  bank.u = u();
  bank.E = settings::run_CE ? E() : g();
  bank.time = time();

  // Convert signed index to a signed surface ID
  if (surface() == SURFACE_NONE) {
    bank.surf_id = SURFACE_NONE;
  } else {
    int surf_id = model::surfaces[surface_index()]->id_;
    bank.surf_id = (surface() > 0) ? surf_id : -surf_id;
  }
}

void Particle::from_source(const SourceSite* src)
{
  // Reset some attributes
  clear();
  surface() = SURFACE_NONE;
  cell_born() = C_NONE;
  material() = C_NONE;
  n_collision() = 0;
  fission() = false;
  zero_flux_derivs();
  lifetime() = 0.0;

  // Copy attributes from source bank site
  type() = src->particle;
  wgt() = src->wgt;
  wgt_last() = src->wgt;
  r() = src->r;
  u() = src->u;
  r_born() = src->r;
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
  parent_nuclide() = src->parent_nuclide;

  // Convert signed surface ID to signed index
  if (src->surf_id != SURFACE_NONE) {
    int index_plus_one = model::surface_map[std::abs(src->surf_id)] + 1;
    surface() = (src->surf_id > 0) ? index_plus_one : -index_plus_one;
  }
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
  if (lowest_coord().cell() == C_NONE) {
    if (!exhaustive_find_cell(*this)) {
      mark_as_lost(
        "Could not find the cell containing particle " + std::to_string(id()));
      return;
    }

    // Set birth cell attribute
    if (cell_born() == C_NONE)
      cell_born() = lowest_coord().cell();

    // Initialize last cells from current cell
    for (int j = 0; j < n_coord(); ++j) {
      cell_last(j) = coord(j).cell();
    }
    n_coord_last() = n_coord();
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
  // Sample a distance to collision
  if (type() == ParticleType::electron || type() == ParticleType::positron) {
    collision_distance() = 0.0;
  } else if (macro_xs().total == 0.0) {
    collision_distance() = INFINITY;
  } else {
    collision_distance() = -std::log(prn(current_seed())) / macro_xs().total;
  }

  // Find the distance to the nearest boundary
  boundary() = distance_to_boundary(*this);

  double speed = this->speed();
  double time_cutoff = settings::time_cutoff[static_cast<int>(type())];
  double distance_cutoff =
    (time_cutoff < INFTY) ? (time_cutoff - time()) * speed : INFTY;

  // Select smaller of the three distances
  double distance =
    std::min({boundary().distance(), collision_distance(), distance_cutoff});

  // Advance particle in space and time
  this->move_distance(distance);
  double dt = distance / speed;
  this->time() += dt;
  this->lifetime() += dt;

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

  // Set particle weight to zero if it hit the time boundary
  if (distance == distance_cutoff) {
    wgt() = 0.0;
  }
}

void Particle::event_cross_surface()
{
  // Saving previous cell data
  for (int j = 0; j < n_coord(); ++j) {
    cell_last(j) = coord(j).cell();
  }
  n_coord_last() = n_coord();

  // Set surface that particle is on and adjust coordinate levels
  surface() = boundary().surface();
  n_coord() = boundary().coord_level();

  if (boundary().lattice_translation()[0] != 0 ||
      boundary().lattice_translation()[1] != 0 ||
      boundary().lattice_translation()[2] != 0) {
    // Particle crosses lattice boundary

    bool verbose = settings::verbosity >= 10 || trace();
    cross_lattice(*this, boundary(), verbose);
    event() = TallyEvent::LATTICE;
  } else {
    // Particle crosses surface
    const auto& surf {model::surfaces[surface_index()].get()};
    // If BC, add particle to surface source before crossing surface
    if (surf->surf_source_ && surf->bc_) {
      add_surf_source_to_bank(*this, *surf);
    }
    this->cross_surface(*surf);
    // If no BC, add particle to surface source after crossing surface
    if (surf->surf_source_ && !surf->bc_) {
      add_surf_source_to_bank(*this, *surf);
    }
    if (settings::weight_window_checkpoint_surface) {
      apply_weight_windows(*this);
    }
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
  surface() = SURFACE_NONE;

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

  if (!model::active_pulse_height_tallies.empty() &&
      type() == ParticleType::photon) {
    pht_collision_energy();
  }

  // Reset banked weight during collision
  n_bank() = 0;
  bank_second_E() = 0.0;
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
    if (coord(j + 1).rotated()) {
      // If next level is rotated, apply rotation matrix
      const auto& m {model::cells[coord(j).cell()]->rotation_};
      const auto& u {coord(j).u()};
      coord(j + 1).u() = u.rotate(m);
    } else {
      // Otherwise, copy this level's direction
      coord(j + 1).u() = coord(j).u();
    }
  }

  // Score flux derivative accumulators for differential tallies.
  if (!model::active_tallies.empty())
    score_collision_derivative(*this);

#ifdef OPENMC_DAGMC_ENABLED
  history().reset();
#endif
}

void Particle::event_revive_from_secondary()
{
  // If particle has too many events, display warning and kill it
  ++n_event();
  if (n_event() == settings::max_particle_events) {
    warning("Particle " + std::to_string(id()) +
            " underwent maximum number of events.");
    wgt() = 0.0;
  }

  // Check for secondary particles if this particle is dead
  if (!alive()) {
    // Write final position for this particle
    if (write_track()) {
      write_particle_track(*this);
    }

    // If no secondary particles, break out of event loop
    if (secondary_bank().empty())
      return;

    from_source(&secondary_bank().back());
    secondary_bank().pop_back();
    n_event() = 0;
    bank_second_E() = 0.0;

    // Subtract secondary particle energy from interim pulse-height results
    if (!model::active_pulse_height_tallies.empty() &&
        this->type() == ParticleType::photon) {
      // Since the birth cell of the particle has not been set we
      // have to determine it before the energy of the secondary particle can be
      // removed from the pulse-height of this cell.
      if (lowest_coord().cell() == C_NONE) {
        bool verbose = settings::verbosity >= 10 || trace();
        if (!exhaustive_find_cell(*this, verbose)) {
          mark_as_lost("Could not find the cell containing particle " +
                       std::to_string(id()));
          return;
        }
        // Set birth cell attribute
        if (cell_born() == C_NONE)
          cell_born() = lowest_coord().cell();

        // Initialize last cells from current cell
        for (int j = 0; j < n_coord(); ++j) {
          cell_last(j) = coord(j).cell();
        }
        n_coord_last() = n_coord();
      }
      pht_secondary_particles();
    }

    // Enter new particle in particle track file
    if (write_track())
      add_particle_track(*this);
  }
}

void Particle::event_death()
{
#ifdef OPENMC_DAGMC_ENABLED
  history().reset();
#endif

  // Finish particle track output.
  if (write_track()) {
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

  if (!model::active_pulse_height_tallies.empty()) {
    score_pulse_height_tally(*this, model::active_pulse_height_tallies);
  }

  // Record the number of progeny created by this particle.
  // This data will be used to efficiently sort the fission bank.
  if (settings::run_mode == RunMode::EIGENVALUE) {
    int64_t offset = id() - 1 - simulation::work_index[mpi::rank];
    simulation::progeny_per_particle[offset] = n_progeny();
  }
}

void Particle::pht_collision_energy()
{
  // Adds the energy particles lose in a collision to the pulse-height

  // determine index of cell in pulse_height_cells
  auto it = std::find(model::pulse_height_cells.begin(),
    model::pulse_height_cells.end(), lowest_coord().cell());

  if (it != model::pulse_height_cells.end()) {
    int index = std::distance(model::pulse_height_cells.begin(), it);
    pht_storage()[index] += E_last() - E();

    // If the energy of the particle is below the cutoff, it will not be sampled
    // so its energy is added to the pulse-height in the cell
    int photon = static_cast<int>(ParticleType::photon);
    if (E() < settings::energy_cutoff[photon]) {
      pht_storage()[index] += E();
    }
  }
}

void Particle::pht_secondary_particles()
{
  // Removes the energy of secondary produced particles from the pulse-height

  // determine index of cell in pulse_height_cells
  auto it = std::find(model::pulse_height_cells.begin(),
    model::pulse_height_cells.end(), cell_born());

  if (it != model::pulse_height_cells.end()) {
    int index = std::distance(model::pulse_height_cells.begin(), it);
    pht_storage()[index] -= E();
  }
}

void Particle::cross_surface(const Surface& surf)
{

  if (settings::verbosity >= 10 || trace()) {
    write_message(1, "    Crossing surface {}", surf.id_);
  }

// if we're crossing a CSG surface, make sure the DAG history is reset
#ifdef OPENMC_DAGMC_ENABLED
  if (surf.geom_type() == GeometryType::CSG)
    history().reset();
#endif

  // Handle any applicable boundary conditions.
  if (surf.bc_ && settings::run_mode != RunMode::PLOTTING) {
    surf.bc_->handle_particle(*this, surf);
    return;
  }

  // ==========================================================================
  // SEARCH NEIGHBOR LISTS FOR NEXT CELL

#ifdef OPENMC_DAGMC_ENABLED
  // in DAGMC, we know what the next cell should be
  if (surf.geom_type() == GeometryType::DAG) {
    int32_t i_cell = next_cell(surface_index(), cell_last(n_coord() - 1),
                       lowest_coord().universe()) -
                     1;
    // save material and temp
    material_last() = material();
    sqrtkT_last() = sqrtkT();
    // set new cell value
    lowest_coord().cell() = i_cell;
    auto& cell = model::cells[i_cell];

    cell_instance() = 0;
    if (cell->distribcell_index_ >= 0)
      cell_instance() = cell_instance_at_level(*this, n_coord() - 1);

    material() = cell->material(cell_instance());
    sqrtkT() = cell->sqrtkT(cell_instance());
    return;
  }
#endif
  int i_surface = std::abs(surface());
  bool verbose = settings::verbosity >= 10 || trace();
  if (surf.is_triso_surface_) {
    if (find_cell_in_virtual_lattice(*this, verbose)) {
      return;
    }
  } else {
    if (neighbor_list_find_cell(*this, verbose)) {
      return;
    }
  }

  // ==========================================================================
  // COULDN'T FIND PARTICLE IN NEIGHBORING CELLS, SEARCH ALL CELLS

  // Remove lower coordinate levels
  n_coord() = 1;
  bool found = exhaustive_find_cell(*this, verbose);

  if (settings::run_mode != RunMode::PLOTTING && (!found)) {
    // If a cell is still not found, there are two possible causes: 1) there is
    // a void in the model, and 2) the particle hit a surface at a tangent. If
    // the particle is really traveling tangent to a surface, if we move it
    // forward a tiny bit it should fix the problem.

    surface() = SURFACE_NONE;
    n_coord() = 1;
    r() += TINY_BIT * u();

    // Couldn't find next cell anywhere! This probably means there is an actual
    // undefined region in the geometry.

    if (!exhaustive_find_cell(*this, verbose)) {
      mark_as_lost("After particle " + std::to_string(id()) +
                   " crossed surface " + std::to_string(surf.id_) +
                   " it could not be located in any cell and it did not leak.");
      return;
    }
  }
}

void Particle::cross_vacuum_bc(const Surface& surf)
{
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

  // Kill the particle
  wgt() = 0.0;

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
  coord(0).cell() = cell_last(0);
  surface() = -surface();

  // If a reflective surface is coincident with a lattice or universe
  // boundary, it is necessary to redetermine the particle's coordinates in
  // the lower universes.
  // (unless we're using a dagmc model, which has exactly one universe)
  n_coord() = 1;
  if (surf.geom_type() != GeometryType::DAG &&
      !neighbor_list_find_cell(*this)) {
    mark_as_lost("Couldn't find particle after reflecting from surface " +
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
  if (settings::max_write_lost_particles < 0 ||
      simulation::n_lost_particles < settings::max_write_lost_particles) {
    write_restart();
  }
  // Increment number of lost particles
  wgt() = 0.0;
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

void Particle::update_neutron_xs(
  int i_nuclide, int i_grid, int i_sab, double sab_frac, double ncrystal_xs)
{
  // Get microscopic cross section cache
  auto& micro = this->neutron_xs(i_nuclide);

  // If the cache doesn't match, recalculate micro xs
  if (this->E() != micro.last_E || this->sqrtkT() != micro.last_sqrtkT ||
      i_sab != micro.index_sab || sab_frac != micro.sab_frac) {
    data::nuclides[i_nuclide]->calculate_xs(i_sab, i_grid, sab_frac, *this);

    // If NCrystal is being used, update micro cross section cache
    if (ncrystal_xs >= 0.0) {
      data::nuclides[i_nuclide]->calculate_elastic_xs(*this);
      ncrystal_update_micro(ncrystal_xs, micro);
    }
  }
}

//==============================================================================
// Non-method functions
//==============================================================================

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

void add_surf_source_to_bank(Particle& p, const Surface& surf)
{
  if (simulation::current_batch <= settings::n_inactive ||
      simulation::surf_source_bank.full()) {
    return;
  }

  // If a cell/cellfrom/cellto parameter is defined
  if (settings::ssw_cell_id != C_NONE) {

    // Retrieve cell index and storage type
    int cell_idx = model::cell_map[settings::ssw_cell_id];

    if (surf.bc_) {
      // Leave if cellto with vacuum boundary condition
      if (surf.bc_->type() == "vacuum" &&
          settings::ssw_cell_type == SSWCellType::To) {
        return;
      }

      // Leave if other boundary condition than vacuum
      if (surf.bc_->type() != "vacuum") {
        return;
      }
    }

    // Check if the cell of interest has been exited
    bool exited = false;
    for (int i = 0; i < p.n_coord_last(); ++i) {
      if (p.cell_last(i) == cell_idx) {
        exited = true;
      }
    }

    // Check if the cell of interest has been entered
    bool entered = false;
    for (int i = 0; i < p.n_coord(); ++i) {
      if (p.coord(i).cell() == cell_idx) {
        entered = true;
      }
    }

    // Vacuum boundary conditions: return if cell is not exited
    if (surf.bc_) {
      if (surf.bc_->type() == "vacuum" && !exited) {
        return;
      }
    } else {

      // If we both enter and exit the cell of interest
      if (entered && exited) {
        return;
      }

      // If we did not enter nor exit the cell of interest
      if (!entered && !exited) {
        return;
      }

      // If cellfrom and the cell before crossing is not the cell of
      // interest
      if (settings::ssw_cell_type == SSWCellType::From && !exited) {
        return;
      }

      // If cellto and the cell after crossing is not the cell of interest
      if (settings::ssw_cell_type == SSWCellType::To && !entered) {
        return;
      }
    }
  }

  SourceSite site;
  site.r = p.r();
  site.u = p.u();
  site.E = p.E();
  site.time = p.time();
  site.wgt = p.wgt();
  site.delayed_group = p.delayed_group();
  site.surf_id = surf.id_;
  site.particle = p.type();
  site.parent_id = p.id();
  site.progeny_id = p.n_progeny();
  int64_t idx = simulation::surf_source_bank.thread_safe_append(site);
}

} // namespace openmc
