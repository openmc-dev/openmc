#include "openmc/ray.h"

#include "openmc/error.h"
#include "openmc/geometry.h"
#include "openmc/material.h"
#include "openmc/mgxs_interface.h"
#include "openmc/settings.h"

namespace openmc {

void Ray::compute_distance()
{
  boundary() = distance_to_boundary(*this);
}

void Ray::trace(double max_distance)
{
  // Clear anything left over from a previous trace on this object.
  reset_trace_state();

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

  // Remaining distance budget
  double max = max_distance;

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

    // Flight through void is real flight: it has to be charged against the
    // distance budget and accumulated, otherwise a ray that starts outside
    // the model reports a traversal distance that excludes the leg up to the
    // model, and max_distance ends up being measured from the entry point
    // rather than from the ray's origin. advance_to_boundary_from_void()
    // leaves boundary().distance() holding the true distance and has already
    // moved the ray that far plus TINY_BIT of padding.
    if (boundary().surface() != SURFACE_NONE &&
        boundary().distance() < INFTY) {
      if (boundary().distance() >= max) {
        // The budget runs out before the ray even reaches the model. Back it
        // up to exactly max, undoing the padding as well.
        move_distance(max - boundary().distance() - TINY_BIT);
        update_distance(max);
        completed_ = true;
        return;
      }
      update_distance(boundary().distance());
      max -= boundary().distance();
    }

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
      warning("Likely infinite loop while ray tracing");
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

    // True geometric distance to the next surface, before any numerical
    // padding is applied. Everything that is physics -- accumulated path
    // length, optical depth, time of flight -- has to be based on this and
    // not on the padded value, or every surface crossing biases the result
    // by TINY_BIT.
    const double surface_distance = boundary().distance();

    // See below comment where call_on_intersection is checked in an
    // if statement for an explanation of this.
    bool call_on_intersection {surface_distance >= 10 * TINY_BIT};

    if (surface_distance >= max) {
      // The ray runs out of budget inside this cell, so no surface is
      // crossed. Only the truncated distance was actually travelled.
      move_distance(max);
      update_distance(max);
      completed_ = true;
      return;
    }

    // DAGMC surfaces expect us to go a little bit further than the advance
    // distance to properly check cell inclusion. The padding is applied to
    // the motion only.
    move_distance(surface_distance + TINY_BIT);

    max -= surface_distance;

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

    // Accumulate before the cell search, while material() still refers to the
    // cell the ray just crossed.
    update_distance(surface_distance);

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
      warning("Likely infinite loop while ray tracing");
      return;
    }
  }
}

void Ray::update_distance(double distance)
{
  // Record how far the ray has traveled
  traversal_distance_ += distance;
}

void ParticleRay::init_physics(
  ParticleType type_, double time_, double E_)
{
  type() = type_;
  time() = time_;
  time_start_ = time_;

  E() = E_;
  E_last() = E_;

  // In multigroup mode the physics is driven by the group index rather than
  // the energy. Particle::from_source() takes the group straight from the
  // source site; here it has to be derived from the energy, otherwise every
  // ray silently uses group 0.
  if (!settings::run_CE) {
    g() = data::mg.get_group_index(E_);
    g_last() = g();
    E() = data::mg.energy_bin_avg_[g()];
    E_last() = E();
  }

  // MacroXS has no default member initializers, and the void branch of
  // update_distance() only writes four of its fields.
  macro_xs() = {};
}

void ParticleRay::mark_as_lost(const char* message)
{
  if (settings::verbosity >= 10) {
    warning(message);
  }
  stop();
}

void ParticleRay::on_intersection() {}

void ParticleRay::update_distance(double distance)
{
  Ray::update_distance(distance);

  time() += distance / speed();

  // Calculate microscopic and macroscopic cross sections
  if (material() != MATERIAL_VOID) {
    if (settings::run_CE) {
      if (material() != material_last() || sqrtkT() != sqrtkT_last() ||
          density_mult() != density_mult_last()) {
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
    macro_xs() = {};
  }

  traversal_mfp_ += macro_xs().total * distance;
}

} // namespace openmc
