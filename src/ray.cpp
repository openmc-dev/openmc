#include "openmc/ray.h"

#include "openmc/error.h"
#include "openmc/geometry.h"
#include "openmc/settings.h"

namespace openmc {

void Ray::compute_distance()
{
  boundary() = distance_to_boundary(*this);
}

void Ray::trace()
{
  // To trace the ray from its origin all the way through the model, we have
  // to proceed in two phases. In the first, the ray may or may not be found
  // inside the model. If the ray is already in the model, phase one can be
  // skipped. Otherwise, the ray has to be advanced to the boundary of the
  // model where all the cells are defined. Importantly, this is assuming that
  // the model is convex, which is a very reasonable assumption for any
  // radiation transport model.
  //
  // After phase one is done, we can starting tracing from cell to cell within
  // the model. This step can use neighbor lists to accelerate the ray tracing.

  bool inside_cell;
  // Check for location if the particle is already known
  if (lowest_coord().cell() == C_NONE) {
    // The geometry position of the particle is either unknown or outside of the
    // edge of the model.
    if (lowest_coord().universe() == C_NONE) {
      // Attempt to initialize the particle. We may have to
      // enter a loop to move it up to the edge of the model.
      inside_cell = exhaustive_find_cell(*this, settings::verbosity >= 10);
    } else {
      // It has been already calculated that the current position is outside of
      // the edge of the model.
      inside_cell = false;
    }
  } else {
    // Availability of the cell means that the particle is located inside the
    // edge.
    inside_cell = true;
  }

  // Advance to the boundary of the model
  while (!inside_cell) {
    advance_to_boundary_from_void();
    inside_cell = exhaustive_find_cell(*this, settings::verbosity >= 10);

    // If true this means no surface was intersected. See cell.cpp and search
    // for numeric_limits to see where we return it.
    if (surface() == std::numeric_limits<int>::max()) {
      warning(fmt::format("Lost a ray, r = {}, u = {}", r(), u()));
      return;
    }

    // Exit this loop and enter into cell-to-cell ray tracing (which uses
    // neighbor lists)
    if (inside_cell)
      break;

    // if there is no intersection with the model, we're done
    if (boundary().surface() == SURFACE_NONE)
      return;

    event_counter_++;
    if (event_counter_ > MAX_INTERSECTIONS) {
      warning("Likely infinite loop in ray traced plot");
      return;
    }
  }

  // Call the specialized logic for this type of ray. This is for the
  // intersection for the first intersection if we had one.
  if (boundary().surface() != SURFACE_NONE) {
    // set the geometry state's surface attribute to be used for
    // surface normal computation
    surface() = boundary().surface();
    on_intersection();
    if (stop_)
      return;
  }

  // reset surface attribute to zero after the first intersection so that it
  // doesn't perturb surface crossing logic from here on out
  surface() = 0;

  // This is the ray tracing loop within the model. It exits after exiting
  // the model, which is equivalent to assuming that the model is convex.
  // It would be nice to factor out the on_intersection at the end of this
  // loop and then do "while (inside_cell)", but we can't guarantee it's
  // on a surface in that case. There might be some other way to set it
  // up that is perhaps a little more elegant, but this is what works just
  // fine.
  while (true) {

    compute_distance();

    // There are no more intersections to process
    // if we hit the edge of the model, so stop
    // the particle in that case. Also, just exit
    // if a negative distance was somehow computed.
    if (boundary().distance() == INFTY || boundary().distance() == INFINITY ||
        boundary().distance() < 0) {
      return;
    }

    // See below comment where call_on_intersection is checked in an
    // if statement for an explanation of this.
    bool call_on_intersection {true};
    if (boundary().distance() < 10 * TINY_BIT) {
      call_on_intersection = false;
    }

    // DAGMC surfaces expect us to go a little bit further than the advance
    // distance to properly check cell inclusion.
    boundary().distance() += TINY_BIT;

    // Advance particle, prepare for next intersection
    for (int lev = 0; lev < n_coord(); ++lev) {
      coord(lev).r() += boundary().distance() * coord(lev).u();
    }
    surface() = boundary().surface();
    // Initialize last cells from the current cell, because the cell() variable
    // does not contain the data for the case of a single-segment ray
    for (int j = 0; j < n_coord(); ++j) {
      cell_last(j) = coord(j).cell();
    }
    n_coord_last() = n_coord();
    n_coord() = boundary().coord_level();
    if (boundary().lattice_translation()[0] != 0 ||
        boundary().lattice_translation()[1] != 0 ||
        boundary().lattice_translation()[2] != 0) {
      cross_lattice(*this, boundary(), settings::verbosity >= 10);
    }

    // Record how far the ray has traveled
    traversal_distance_ += boundary().distance();
    inside_cell = neighbor_list_find_cell(*this, settings::verbosity >= 10);

    // Call the specialized logic for this type of ray. Note that we do not
    // call this if the advance distance is very small. Unfortunately, it seems
    // darn near impossible to get the particle advanced to the model boundary
    // and through it without sometimes accidentally calling on_intersection
    // twice. This incorrectly shades the region as occluded when it might not
    // actually be. By screening out intersection distances smaller than a
    // threshold 10x larger than the scoot distance used to advance up to the
    // model boundary, we can avoid that situation.
    if (call_on_intersection) {
      on_intersection();
      if (stop_)
        return;
    }

    if (!inside_cell)
      return;

    event_counter_++;
    if (event_counter_ > MAX_INTERSECTIONS) {
      warning("Likely infinite loop in ray traced plot");
      return;
    }
  }
}

} // namespace openmc
