#include "openmc/ray.h"

#include "openmc/error.h"
#include "openmc/geometry.h"
#include "openmc/material.h"
#include "openmc/mgxs_interface.h"
#include "openmc/settings.h"

namespace openmc {

namespace {
// Max intersections before we assume ray tracing is caught in an infinite loop
constexpr int MAX_INTERSECTIONS = 1000000;
} // namespace

template<typename RayT>
void trace_ray(RayT& ray, double max_distance)
{
  // Clear anything left over from a previous trace on this object.
  ray.reset_trace_state();

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
  if (ray.lowest_coord().cell() == C_NONE) {
    // The geometry position of the particle is either unknown or outside of the
    // edge of the model.
    if (ray.lowest_coord().universe() == C_NONE) {
      // Attempt to initialize the particle. We may have to
      // enter a loop to move it up to the edge of the model.
      inside_cell = exhaustive_find_cell(ray, settings::verbosity >= 10);
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
    ray.advance_to_boundary_from_void();

    // Flight through void is real flight: it has to be charged against the
    // distance budget, otherwise max_distance ends up being measured from the
    // point where the ray entered the model rather than from its origin. It
    // is accumulated into total_distance_ only -- traversal_distance_ stays
    // zero until the model is reached.
    if (ray.boundary().surface() != SURFACE_NONE &&
        ray.boundary().distance() < INFTY) {
      // advance_to_boundary_from_void() has already moved the ray to the
      // boundary plus the TINY_BIT of padding, so that is the distance to
      // account for.
      double advance = ray.boundary().distance() + TINY_BIT;

      if (advance >= max) {
        // The budget runs out before the ray even reaches the model. Back it
        // up so the net movement is exactly max.
        ray.move_distance(max - advance);
        ray.update_distance(max);
        ray.completed_ = true;
        return;
      }
      ray.update_distance(advance);
      max -= advance;
    }

    inside_cell = exhaustive_find_cell(ray, settings::verbosity >= 10);

    // If true this means no surface was intersected. See cell.cpp and search
    // for numeric_limits to see where we return it.
    if (ray.surface() == std::numeric_limits<int>::max()) {
      warning(fmt::format("Lost a ray, r = {}, u = {}", ray.r(), ray.u()));
      return;
    }

    // Exit this loop and enter into cell-to-cell ray tracing (which uses
    // neighbor lists)
    if (inside_cell)
      break;

    // if there is no intersection with the model, we're done
    if (ray.boundary().surface() == SURFACE_NONE)
      return;

    ray.event_counter_++;
    if (ray.event_counter_ > MAX_INTERSECTIONS) {
      warning("Likely infinite loop while ray tracing");
      return;
    }
  }

  // From here on the ray is inside the model, so its flight counts toward
  // traversal_distance_ as well as total_distance_.
  ray.in_model_ = true;

  // Call the specialized logic for this type of ray. This is for the
  // intersection for the first intersection if we had one.
  if (ray.boundary().surface() != SURFACE_NONE) {
    // set the geometry state's surface attribute to be used for
    // surface normal computation
    ray.surface() = ray.boundary().surface();
    ray.on_intersection();
    if (ray.stop_)
      return;
  }

  // reset surface attribute to zero after the first intersection so that it
  // doesn't perturb surface crossing logic from here on out
  ray.surface() = 0;

  // This is the ray tracing loop within the model. It exits after exiting
  // the model, which is equivalent to assuming that the model is convex.
  // It would be nice to factor out the on_intersection at the end of this
  // loop and then do "while (inside_cell)", but we can't guarantee it's
  // on a surface in that case. There might be some other way to set it
  // up that is perhaps a little more elegant, but this is what works just
  // fine.
  while (true) {

    ray.boundary() = distance_to_boundary(ray);

    // There are no more intersections to process
    // if we hit the edge of the model, so stop
    // the particle in that case. Also, just exit
    // if a negative distance was somehow computed.
    if (ray.boundary().distance() == INFTY ||
        ray.boundary().distance() == INFINITY ||
        ray.boundary().distance() < 0) {
      return;
    }

    // Distance from the ray's current position to the next surface.
    const double surface_distance = ray.boundary().distance();

    // See below comment where call_on_intersection is checked in an
    // if statement for an explanation of this.
    bool call_on_intersection {surface_distance >= 10 * TINY_BIT};

    // DAGMC surfaces expect us to go a little bit further than the advance
    // distance to properly check cell inclusion.
    double advance = surface_distance + TINY_BIT;

    if (advance >= max) {
      // The ray runs out of budget inside this cell, so no surface is
      // crossed. Only the truncated distance was actually travelled.
      ray.move_distance(max);
      ray.update_distance(max);
      ray.completed_ = true;
      return;
    }

    ray.move_distance(advance);

    max -= advance;

    ray.surface() = ray.boundary().surface();
    // Initialize last cells from the current cell, because the cell() variable
    // does not contain the data for the case of a single-segment ray
    for (int j = 0; j < ray.n_coord(); ++j) {
      ray.cell_last(j) = ray.coord(j).cell();
    }
    ray.n_coord_last() = ray.n_coord();
    ray.n_coord() = ray.boundary().coord_level();
    if (ray.boundary().lattice_translation()[0] != 0 ||
        ray.boundary().lattice_translation()[1] != 0 ||
        ray.boundary().lattice_translation()[2] != 0) {
      cross_lattice(ray, ray.boundary(), settings::verbosity >= 10);
    }

    // Accumulate before the cell search, while material() still refers to the
    // cell the ray just crossed.
    ray.update_distance(advance);

    inside_cell = neighbor_list_find_cell(ray, settings::verbosity >= 10);

    // Call the specialized logic for this type of ray. Note that we do not
    // call this if the advance distance is very small. Unfortunately, it seems
    // darn near impossible to get the particle advanced to the model boundary
    // and through it without sometimes accidentally calling on_intersection
    // twice. This incorrectly shades the region as occluded when it might not
    // actually be. By screening out intersection distances smaller than a
    // threshold 10x larger than the scoot distance used to advance up to the
    // model boundary, we can avoid that situation.
    if (call_on_intersection) {
      ray.on_intersection();
      if (ray.stop_)
        return;
    }

    if (!inside_cell)
      return;

    ray.event_counter_++;
    if (ray.event_counter_ > MAX_INTERSECTIONS) {
      warning("Likely infinite loop while ray tracing");
      return;
    }
  }
}

void RayState::accumulate_distance(double distance)
{
  total_distance_ += distance;
  if (in_model_) {
    traversal_distance_ += distance;
  }
}

void Ray::trace(double max_distance)
{
  trace_ray(*this, max_distance);
}

void ParticleRay::trace(double max_distance)
{
  trace_ray(*this, max_distance);
}

void ParticleRay::init_physics(ParticleType type_, double time_, double E_)
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

void ParticleRay::update_distance(double distance)
{
  accumulate_distance(distance);

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

// Explicit instantiations: the two kinds of ray that share the tracing loop
template void trace_ray<Ray>(Ray&, double);
template void trace_ray<ParticleRay>(ParticleRay&, double);

} // namespace openmc
