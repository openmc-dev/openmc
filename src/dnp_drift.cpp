#include "openmc/dnp_drift.h"
#include "openmc/error.h"
#include "openmc/field.h"
#include "openmc/particle_data.h"
#include "openmc/position.h"
#include "openmc/simulation.h"

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
    fatal_error(
      "Time t at point B must be greater than tine ta at point A.");
  }

  double ab = (pb - pa).norm();
  double ac = (pc - pa).norm();
  double cb = (pb - pc).norm();

  if (ab - (ac + cb) > DNP_DRIFT_DISTANCE_MIN) {
    fatal_error(
      "Point C must be located between point A and point B.");
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
  cell_n = simulation::velocity_field.get_mesh_bin(y_n);
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
          action = Actions::PLACE_AT_INLET;
        } else {
          action = Actions::BLOCK_AT_OUTLET;
        }
        break;
      case BCType::INLET:
        action = Actions::PLACE_AT_INLET;
        break;
      case BCType::WALL:
        action = Actions::BLOCK_AT_LOCATION;
        break;
      default:
        fatal_error("Not implemented for this boundary type!");
        break;
      }

      // Perform action
      switch (action) {
      case Actions::PLACE_AT_INLET:
        if (decay_time > t_n + settings::dnp_drift_external_travel_time) {
          // The DNP has time to reenter the modeled part of the system
          simulation::velocity_field.randomly_place_on_inlet(y_n, cell_n, seed);
          t_n += settings::dnp_drift_external_travel_time;
          cell_n_minus_1 = -1;
          y_n_minus_1 = Position();
        } else {
          // The DNP decays outside the modeled part of the system
          t_before_decay = 0.;
          return false;
        }
        break;
      case Actions::BLOCK_AT_OUTLET:
        t_before_decay = 0.;
        return false; // TODO - return true normally
        break;
      case Actions::BLOCK_AT_LOCATION:
        t_before_decay = decay_time - t_n;
        return false; // TODO - return true normally
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

} // namespace openmc
