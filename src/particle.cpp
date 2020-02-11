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
Particle::create_secondary(Direction u, double E, Type type)
{
  secondary_bank_.emplace_back();

  auto& bank {secondary_bank_.back()};
  bank.particle = type;
  bank.wgt = wgt_;
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
    if (!find_cell(this, false)) {
      this->mark_as_lost("Could not find the cell containing particle "
        + std::to_string(id_));
      return;
    }

    // Set birth cell attribute
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
  boundary_ = distance_to_boundary(this);

  // Sample a distance to collision
  if (type_ == Particle::Type::electron ||
      type_ == Particle::Type::positron) {
    collision_distance_ = 0.0;
  } else if (macro_xs_.total == 0.0) {
    collision_distance_ = INFINITY;
  } else {
    collision_distance_ = -std::log(prn(this->current_seed())) / macro_xs_.total;
    // cvmt: continuous varying materials tracking  
    if(model::materials[material_]->continuous_num_density_) 
	      collision_distance_ = sampling_cvmt(this, boundary_.distance);
    //
  }

  // Select smaller of the two distances
  double distance = std::min(boundary_.distance, collision_distance_);

  // Advance particle
  for (int j = 0; j < n_coord_; ++j) {
    coord_[j].r += distance * coord_[j].u;
  }

  // Score track-length tallies
  if (!model::active_tracklength_tallies.empty()) {
    score_tracklength_tally(this, distance);
  }

  // Score track-length estimate of k-eff
  if (settings::run_mode == RunMode::EIGENVALUE &&
      type_ == Particle::Type::neutron) {
    keff_tally_tracklength_ += wgt_ * distance * macro_xs_.nu_fission;
  }

  // Score flux derivative accumulators for differential tallies.
  if (!model::active_tallies.empty()) {
    score_track_derivative(this, distance);
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
    cross_lattice(this, boundary_);
    event_ = TallyEvent::LATTICE;
  } else {
    // Particle crosses surface
    this->cross_surface();
    event_ = TallyEvent::SURFACE;
  }
  // Score cell to cell partial currents
  if (!model::active_surface_tallies.empty()) {
    score_surface_tally(this, model::active_surface_tallies);
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
  if (!model::active_tallies.empty()) score_collision_derivative(this);
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
    write_message("    Crossing surface " + std::to_string(surf->id_));
  }

  if (surf->bc_ == Surface::BoundaryType::VACUUM && (settings::run_mode != RunMode::PLOTTING)) {
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
    keff_tally_leakage_ += wgt_;

    // Display message
    if (settings::verbosity >= 10 || trace_) {
      write_message("    Leaked out of surface " + std::to_string(surf->id_));
    }
    return;

  } else if ((surf->bc_ == Surface::BoundaryType::REFLECT ||
              surf->bc_ == Surface::BoundaryType::WHITE)
              && (settings::run_mode != RunMode::PLOTTING)) {
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

    Direction u = (surf->bc_ == Surface::BoundaryType::REFLECT) ?
      surf->reflect(this->r(), this->u(), this) :
      surf->diffuse_reflect(this->r(), this->u(), this->current_seed());

    // Make sure new particle direction is normalized
    this->u() = u / u.norm();

    // Reassign particle's cell and surface
    coord_[0].cell = cell_last_[n_coord_last_ - 1];
    surface_ = -surface_;

    // If a reflective surface is coincident with a lattice or universe
    // boundary, it is necessary to redetermine the particle's coordinates in
    // the lower universes.
    // (unless we're using a dagmc model, which has exactly one universe)
    if (!settings::dagmc) {
      n_coord_ = 1;
      if (!find_cell(this, true)) {
        this->mark_as_lost("Couldn't find particle after reflecting from surface "
                           + std::to_string(surf->id_) + ".");
        return;
      }
    }

    // Set previous coordinate going slightly past surface crossing
    r_last_current_ = this->r() + TINY_BIT*this->u();

    // Diagnostic message
    if (settings::verbosity >= 10 || trace_) {
      write_message("    Reflected from surface " + std::to_string(surf->id_));
    }
    return;

  } else if (surf->bc_ == Surface::BoundaryType::PERIODIC && settings::run_mode != RunMode::PLOTTING) {
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
    if (settings::verbosity >= 10 || trace_) {
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

  if (settings::run_mode != RunMode::PLOTTING && (!found)) {
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
  if (settings::run_mode == RunMode::PARTICLE) return;

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
      case RunMode::FIXED_SOURCE:
        write_dataset(file_id, "run_mode", "fixed source");
        break;
      case RunMode::EIGENVALUE:
        write_dataset(file_id, "run_mode", "eigenvalue");
        break;
      case RunMode::PARTICLE:
        write_dataset(file_id, "run_mode", "particle restart");
        break;
    }
    write_dataset(file_id, "id", id_);
    write_dataset(file_id, "type", static_cast<int>(type_));

    int64_t i = current_work_;
    write_dataset(file_id, "weight", simulation::source_bank[i-1].wgt);
    write_dataset(file_id, "energy", simulation::source_bank[i-1].E);
    write_dataset(file_id, "xyz", simulation::source_bank[i-1].r);
    write_dataset(file_id, "uvw", simulation::source_bank[i-1].u);

    // Close file
    file_close(file_id);
  } // #pragma omp critical
}

//==============================================================================
// cvmt implementation 
//==============================================================================
double 
Particle::sampling_cvmt(Particle* p, double d_boundary)
{
  double PNC, tau_hat;            
  double d_collision {0.0};
  double optical_depth {0.0};
  std::vector<double> xs_t; 
  Position xyz_orig;  
  //
  xs_t.resize(settings::num_intervals + 1);
  
  simpsons_path_integration(optical_depth, d_boundary, xs_t, false, 0);
  
  PNC = std::exp(-optical_depth);
  
  if(prn(this->current_seed()) <= PNC) { // no collision
    d_collision = INFINITY;
  } else { // collision   
    // sample optical depth 
    tau_hat = -std::log(1.0 - (1.0 - PNC) * prn(this->current_seed()));
    // get the flight distance for sampled optical depth
    estimate_flight_distance(xs_t, d_boundary, tau_hat, d_collision);
    // move particle to the point of collision so we can make sure that
    xyz_orig = p->coord_[p->n_coord_].r;
    // cross sections are updated for later tallying of keff and user tallies
    move_particle_coord(p->coord_[p->n_coord_], d_collision);
    model::materials[material_]->calculate_xs(*p);
    // move particle back
    p->coord_[p->n_coord_].r = xyz_orig;
  }
  return d_collision;
}

void 
Particle::move_particle_coord(LocalCoord& coord, double ds) {coord.r += ds * coord.u;}

void 
Particle::copy_data(LocalCoord& to, LocalCoord from)
{
  to.r         = from.r;
  to.u         = from.u;
  to.cell      = from.cell;
  to.universe  = from.universe;
  to.lattice   = from.lattice;
  to.lattice_x = from.lattice_x;
  to.lattice_y = from.lattice_y;
  to.lattice_z = from.lattice_z;
  to.rotated   = from.rotated;
}

void 
Particle::simpsons_path_integration(double& optical_depth, double distance, std::vector<double>& xs_t, bool dbg_file, int it_num)
{
  double ds, new_temp;
  LocalCoord coord;
  
  copy_data(coord, this->coord_[this->n_coord_ - 1]);
  
  ds = distance / (double)(settings::num_intervals);

  for (int i = 0; i < settings::num_intervals + 1; i++){ 
    //! Recalculate the cross section
    model::materials[material_]->calculate_xs(*this);
    //! Save the total cross section
    xs_t[i] = macro_xs_.total;   
    //! Move particle along path
    move_particle_coord(this->coord_[this->n_coord_ - 1], ds);
  }
  //
  copy_data(this->coord_[this -> n_coord_ - 1], coord);
  //
  optical_depth = 0.0; 
  for(int i = 0; i < settings::num_intervals - 1; i = i + 2){
     //! Add to the integral
     optical_depth += 2.0 * ds / 6.0 * ( xs_t[i] + 4.0 * xs_t[i+1] + xs_t[i+2]);
  }
  return; 
}

void 
Particle::estimate_flight_distance(std::vector<double> xs_t, double distance, double tau_hat, double& s)
{
  int index {1};
  double a, b, c, d, m, dds;
  double ds, delta_tau_hat, optical_depth;
  bool tau_overrun {false};
  //! Calculate differential path length
  ds = distance / (double)(settings::num_intervals);
  //! Set variables for loop
  optical_depth = 0.0;
  //! Find the index of the cell that goes over the sampled path integration
  while (optical_depth <= tau_hat) {
    if (index > (settings::num_intervals - 1)) {
      tau_overrun = true;
      optical_depth = tau_hat;
    } else {
      optical_depth = optical_depth + \
                      2.0 * ds / 6.0 * ( xs_t[index - 1] + \
                      4.0 * xs_t[index] + xs_t[index + 1]);
      index += 2;
    }
  }
  //! Subtract off the last delta optical depth that we added
  index -= 2;
  optical_depth = optical_depth - \
                  2.0 * ds / 6.0 * (xs_t[index - 1] + \
                  4.0 * xs_t[index] + xs_t[index + 1]);
  //! Calculate the delta in the optical depth
  delta_tau_hat = tau_hat - optical_depth;
  //
  if (tau_overrun) {
    warning("The search for the optical depth bin failed.");
    s = distance;
  } else {
    //! Determine the coefficients for a second order expansion in that bin
    //! sigma_t = ax^2 + bx + c
    //! Note that we are shifting the portion of the curve so that the first
    //! point lies at zero.  This makes integration much easier.
    c = xs_t[index - 1];
    b = (2.0 * xs_t[index] - 0.5 * xs_t[index + 1] - 1.5 * c) / ds;
    a = (xs_t[index + 1] - c - 2.0 * ds * b) / (4.0 * ds * ds);
    //! Define quantities for solution of the cubic equation
    b /= 2.0;
    a /= 3.0;
    d = -delta_tau_hat;
    //! If the polynomial is not truly cubic, the soluiton will blow up
    //! Choose between integrated polynomial of cubic,quadratic, and linear
    if (std::abs(a) > 1E-10) {
      get_cubic_root(a, b, c, d, 0.0, 2.0 * ds, dds);
      s = ds * (index - 1) + dds;
    } else if (std::abs(b) > 1E-10) {
      //! Get the best guest for root based on a linear fit
      m = (xs_t[index + 1] - xs_t[index - 1]) / (2.0 * ds);
      get_quadratic_root(0.5 * m, c, d, 0.0, 2.0 * ds, dds);
      s = ds * (index-1) + dds;
    } else {
      //! This is the case where the cross section is constant.
      s = ds * (index - 1) - d / c;
      //! Put some checks in for now
      if (s < (ds * (index - 1)) || s > (ds * (index + 1))) {
        fatal_error("s for the constant cross section case lies out of bounds.");
      }
    }
  }
  if (s > distance || s < 0.0) {
    fatal_error("The predicted distance is greater than distance to boundary or less than zero.");
  }
  return;
}

double
Particle::sign(double aa, double bb)
{
  if(aa * bb > 0.0) 
      return aa;
  else 
      return -aa;
}

void 
Particle::get_quadratic_root(double a, double b, double c, double lower_b, double upper_b, double& root)
{
  double C_1, C_2;
  double e1, e2, e3, e4;
  double q, min_error, inner_sqrt;
  //
  inner_sqrt = b * b - 4.0 * a * c;
  
  if (inner_sqrt < 0.0) {
    fatal_error("Non real roots for quadratic equation.");
  } else if (inner_sqrt >= 0) {  
    q = -1.0 / 2.0 * (b + sign(1.0,b) * std::sqrt(inner_sqrt));
    C_1 = q / a;
    C_2 = c / q ; 
    if (C_1 >= lower_b && C_1 <= upper_b) {
      root = C_1;
    } else if (C_2 >= lower_b && C_2 <= upper_b) {
      root = C_2;
    } else {
      e1 = std::abs(C_1 - lower_b);
      e2 = std::abs(C_1 - upper_b);
      e3 = std::abs(C_2 - lower_b);
      e4 = std::abs(C_2 - upper_b);
      min_error = std::min(std::min(std::min(e1, e2), e3), e4);
      if (std::abs(e1) == min_error || std::abs(e3) == min_error) {
        root = lower_b;
      } else if (std::abs(e2) == min_error || std::abs(e4) == min_error) {
        root = upper_b;
      } else {
        fatal_error("Invalid case for handling quadratic bounds overload.");
      }  
      warning("Neither quadratic roots satisifed the criteria.");
      warning("Quadratic roots not within 1E-5 of bounds.");
    }
  }
  return;
}

void 
Particle::get_cubic_root(double a, double b, double c, double d, \
		         double lower_b, double upper_b, double& root)
{
  double aa, bb, cc;  // Coefficients in x^3 + aa * x^2 + bb * x + cc = 0 
  double Q, R, dd, ee, theta;
  // Calculate coefficients as they appear in the numerical recipes book
  aa = b / a;
  bb = c / a;
  cc = d / a;
  Q = (aa * aa - 3.0 * bb) / 9.0;
  R = (2.0 * aa * aa * aa - 9.0 * aa * bb + 27.0 * cc) / 54.0;
  if (R*R < Q*Q*Q) {
    theta = std::acos(R / std::sqrt(Q*Q*Q));
    root = -2.0 * std::sqrt(Q) * std::cos(theta/3.0) - aa / 3.0;
    if (root < lower_b || root > upper_b) {
        root = -2.0 * std::sqrt(Q) * std::cos( (theta + 2.0 * PI) / 3.0) - aa / 3.0;
        if (root < lower_b || root > upper_b) {
           root = -2.0 * std::sqrt(Q) * std::cos((theta - 2.0 * PI) / 3.0) - aa / 3.0;
           if (root < lower_b || root > upper_b) {
              if (root < lower_b) {
                 root = lower_b;
              } else {
                 root = upper_b;
              }
              warning("Acceptable cubic root not found.");
           }
        }
    }
  } else {
    if (R >= 0.0) {
        dd = -std::pow((std::abs(R) + std::sqrt(std::abs(R * R - Q * Q * Q))), 1.0 / 3.0);
    } else {
        dd = std::pow(std::abs(R) + std::sqrt(std::abs(R * R - Q * Q * Q)), 1.0 / 3.0);
    } 
    if (dd == 0.0) {
        ee = 0.0;
    } else {
        ee = Q / dd;
    }
    root = (dd + ee) - aa / 3.0;
    if (root < lower_b || root > upper_b) {     
      if (root < lower_b) {
           root = lower_b;
      } else {
           root = upper_b;
      }
      warning("Acceptable cubic root not found.");
    }
  }
  return; 
}

} // namespace openmc
