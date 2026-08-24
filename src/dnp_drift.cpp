#include "openmc/dnp_drift.h"
#include "openmc/cell.h"
#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/field.h"
#include "openmc/geometry.h"
#include "openmc/particle_data.h"
#include "openmc/position.h"
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
  cell_n = simulation::velocity_field.get_bin(y_n);
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
    cell_n = simulation::velocity_field.get_next_bin(
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
        simulation::velocity_field.randomly_place_on_inlet(y_n, cell_n, seed);
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

bool reconcile_precursor_drift(SourceSite& site)
{
  Particle p;
  p.r() = site.r;
  p.u() = site.u;

  // Is the DNP inside the model?
  bool found = exhaustive_find_cell(p);

  // If outside, it might be because we are on a boundary with an outward
  // direction
  if (!found) {

    // Nudge the particle backward
    p.r() = p.r() - p.u() * TINY_BIT;
    found = exhaustive_find_cell(p);

    // If not found here, it is certainly lost
    if (!found)
      fatal_error("DNP is certainly lost");

    // Go back to the previous position
    p.r() = site.r;
  }

  // Apply boundary condition if the site is on a surface
  for (auto& surf_id : model::cells[p.lowest_coord().cell()]->surfaces()) {
    const Surface& surf = *model::surfaces[std::abs(surf_id) - 1].get();
    double eval = surf.evaluate(p.r());
    if (std::abs(eval) < FP_COINCIDENT) {

      Direction normal = surf.normal(p.r());
      double dot = site.u.dot(normal);
      bool going_outward = dot > 0.0;

      // TODO: manage albedo coefficient
      // TODO: add more boundary types

      // Transmission
      if (surf.bc_->type() == "transmission")
        return true;

      // Vacuum
      if (surf.bc_->type() == "vacuum")
        return (!going_outward);

      // Reflective
      if (surf.bc_->type() == "reflective") {
        if (going_outward) {
          // Apply reflection
          site.u = surf.reflect(p.r(), p.u(), &p);
          site.u /= site.u.norm();
        }
        return true;
      }

      // Other
      fatal_error("Not implemented!");
    }
  }

  return found;
}

} // namespace openmc
