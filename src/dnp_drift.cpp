#include "openmc/dnp_drift.h"
#include "openmc/cell.h"
#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/field.h"
#include "openmc/geometry.h"
#include "openmc/particle_data.h"
#include "openmc/position.h"
#include "openmc/random_lcg.h"
#include "openmc/simulation.h"
#include "openmc/surface.h"

namespace openmc {

void _adjust_position(Position& y_n, const Position& y_n_minus_1, double t,
  double dt, double decay_time)
{
  if (dt <= 0.0) {
    fatal_error("The time step dt should be greater than 0.");
  }
  if (decay_time > t) {
    fatal_error("Decay time should be lower or equal to the total time t.");
  }

  double excess_time = t - decay_time;

  if (excess_time > dt) {
    fatal_error("The excess time cannot be greater than the time step dt.");
  }

  double ab = (y_n - y_n_minus_1).norm();

  if (ab <= 0.0) {
    fatal_error(
      "The distance between the two points cannot be lower or equal to 0.");
  }

  double ac = (dt - excess_time) * ab / dt;

  y_n = y_n_minus_1 + (y_n - y_n_minus_1) * ac / ab;
}

void _adjust_time(double& t, double ta, const Position& pa, const Position& pb,
  const Position& pc)
{
  if (t <= ta) {
    fatal_error("Time t at point B must be greater than tine ta at point A.");
  }

  double ab = (pb - pa).norm();
  double ac = (pc - pa).norm();
  double cb = (pb - pc).norm();

  if (ab - (ac + cb) > DNP_DRIFT_DISTANCE_MIN) {
    fatal_error("Point C must be located between point A and point B.");
  }

  t -= (t - ta) * cb / ab;
}

bool transport_dnp(SourceSite& site, double decay_time, uint64_t* seed)
{
  // Initialization
  double t_n = 0.0;                   // Current (n) integration time
  double t_n_minus_1;                 // Previous (n-1) integration time
  Position& y_n = site.r;             // Current (n) location of the site
  Position y_n_minus_1;               // Previous (n-1) location of the site
  int cell_n;                         // Current (n) cell
  int cell_n_minus_1;                 // Previous (n-1) cell
  int n_iteration = 0;                // Total number of performed iterations
  double t_before_decay = decay_time; // Time remaining before decay occurs
  BCType crossed_boundary =
    BCType::NONE;        // Boundary type of the last crossed surface
  Position intersection; // Intersection with the mesh boundary

  // Localize initial position
  cell_n = simulation::velocity_field->get_bin(y_n);
  if (cell_n < 0) {
    //  Particle can be close to the wall of the mesh so it can be stopped at
    //  the initial location because we generally put no slip boundary
    //  conditions on walls which means that the particle would not have been
    //  able to move even if detected inside the mesh
    t_before_decay = 0.;
    return true;
  }

  // Transport
  while (t_n < decay_time) {

    // Increment number of iterations
    n_iteration++;

    // Verify that we have not reached the maximum number of iterations
    if (n_iteration > DNP_DRIFT_TRANSPORT_MAX_ITER) {
      fatal_error("Maximum number of iterations reached during DNP transport!");
    }

    // Save previous states before integration
    t_n_minus_1 = t_n;
    y_n_minus_1 = y_n;
    cell_n_minus_1 = cell_n;

    // Integration
    simulation::streamline_integrator->next_step(
      t_n, y_n, cell_n, simulation::velocity_field);

    // Find the next cell
    cell_n = simulation::velocity_field->get_next_bin(
      y_n_minus_1, y_n, cell_n, crossed_boundary, intersection);

    // If the distance between two consecutive points is low, we block the point
    // in position
    if ((y_n - y_n_minus_1).norm() < DNP_DRIFT_DISTANCE_MIN) {
      t_before_decay = 0.;
      return true;
    }

    // If the particle left the mesh
    if (cell_n < 0) {

      // Update time and position
      _adjust_time(t_n, t_n_minus_1, y_n_minus_1, y_n, intersection);
      y_n = intersection;

      // If decay time occurs before intersection, stop when decay time is
      // reached
      if (t_n > decay_time) {
        _adjust_position(y_n, y_n_minus_1, t_n, t_n - t_n_minus_1, decay_time);
        t_n = decay_time;
        t_before_decay = 0.;
        return true;
      }

      // Determine the action to perform based on the crossed boundary
      Actions action;

      switch (crossed_boundary) {
      case BCType::OUTLET:
        if (settings::dnp_drift_recycling_on) {
          if (decay_time > t_n + settings::dnp_drift_external_travel_time) {
            action = Actions::REENTER;
          } else {
            action = Actions::ESCAPE;
          }
        } else {
          action = Actions::ESCAPE;
        }
        break;
      case BCType::INLET:
        action = Actions::REENTER;
        break;
      case BCType::WALL:
        action = Actions::DECAY_IN_PLACE;
        break;
      default:
        fatal_error("Not implemented for this boundary type!");
        break;
      }

      // Perform action
      switch (action) {
      case Actions::REENTER:
        simulation::velocity_field->randomly_place_on_inlet(y_n, cell_n, seed);
        t_n += settings::dnp_drift_external_travel_time;
        cell_n_minus_1 = -1;
        y_n_minus_1 = Position();
        break;
      case Actions::ESCAPE:
        t_before_decay = decay_time - t_n;
        return false;
        break;
      case Actions::DECAY_IN_PLACE:
        t_before_decay = 0.;
        return true;
        break;
      default:
        fatal_error("Unrecognized action in DNP transport!");
        break;
      }
    }
  }

  // Adjust time and position to fit exactly
  if (t_n > decay_time) {
    _adjust_position(y_n, y_n_minus_1, t_n,
      simulation::streamline_integrator->dt(), decay_time);
    t_n = decay_time;
  }

  t_before_decay = 0;
  return true;
}

bool reconcile_precursor_drift(SourceSite& site, int64_t particle_id)
{
  // Set up temporary particle at the precursor's position and direction
  Particle p;
  int64_t particle_seed = compute_transport_seed(particle_id);
  init_particle_seeds(particle_seed, p.seeds());
  p.stream() = STREAM_TRACKING;
  p.r() = site.r;
  p.u() = site.u;

  // Case 1: the particle is seen inside the geometry.
  if (exhaustive_find_cell(p)) {
    return true;
  }

  // Case 2: the particle is seen outside the geometry. This means that the
  // particle is genuinely outside the geometry, on a boundary but heading
  // outward, or traveling exactly along a boundary plane.

  // Nudge the particle slightly backward to check if it reenters.
  p.r() = p.r() - p.u() * TINY_BIT;

  if (!exhaustive_find_cell(p)) {
    // The particle is not recoverable: it is either in an undefined region or
    // traveling exactly along a surface plane.
    warning("DNP lost: unable to locate cell after backward nudge. If this "
            "warning is occurring frequently, this might indicate spatial "
            "inconsistency between the mesh used for DNP transport and the "
            "OpenMC model. Please check both the mesh and the OpenMC model.");
    return false;
  }

  // Compute the distance from the nudged position to the nearest boundary.
  BoundaryInfo boundary = distance_to_boundary(p);

  // If the boundary is closer than the nudge distance, the particle was
  // genuinely outside before the nudge and not on a boundary.
  if (boundary.distance() < TINY_BIT - FP_COINCIDENT) {
    return false;
  }

  // The particle is on a boundary heading outward.
  // Reset particle to its original state and apply boundary conditions.
  p.r() = site.r;
  p.u() = site.u;

  const auto& surf = *model::surfaces[std::abs(boundary.surface()) - 1];
  const auto* bc = surf.bc_.get();

  if (!bc) {
    fatal_error(
      "Boundary condition not found after crossing an external boundary!");
  }

  const std::string& bc_type = bc->type();

  // Transmission
  if (bc_type == "transmission") {
    fatal_error("Potentially crossing an external transmission boundary!");
  }

  // Vacuum
  if (bc_type == "vacuum") {
    return false;
  }

  // Remaining conditions: reflective, white, or periodic
  // Apply transformation and albedo
  bc->transform(p, surf, site.r, site.u, site.surf_id);
  if (bc->has_albedo()) {
    site.wgt *= bc->albedo();
  }

  return true;
}

} // namespace openmc
