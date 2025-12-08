#include "openmc/boundary_condition.h"

#include <exception>

#include <fmt/core.h>

#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/random_ray/random_ray.h"
#include "openmc/surface.h"

namespace openmc {

//==============================================================================
// VacuumBC implementation
//==============================================================================

void VacuumBC::handle_particle(Particle& p, const Surface& surf) const
{
  // Random ray and Monte Carlo need different treatments at vacuum BCs
  if (settings::solver_type == SolverType::RANDOM_RAY) {
    // Reflect ray off of the surface
    ReflectiveBC().handle_particle(p, surf);

    // Set ray's angular flux spectrum to vacuum conditions (zero)
    RandomRay* r = static_cast<RandomRay*>(&p);
    std::fill(r->angular_flux_.begin(), r->angular_flux_.end(), 0.0);

  } else {
    p.cross_vacuum_bc(surf);
  }
}

//==============================================================================
// ReflectiveBC implementation
//==============================================================================

void ReflectiveBC::handle_particle(Particle& p, const Surface& surf) const
{
  Direction u = surf.reflect(p.r(), p.u(), &p);
  u /= u.norm();

  // Handle the effects of the surface albedo on the particle's weight.
  BoundaryCondition::handle_albedo(p, surf);

  p.cross_reflective_bc(surf, u);
}

//==============================================================================
// WhiteBC implementation
//==============================================================================

void WhiteBC::handle_particle(Particle& p, const Surface& surf) const
{
  Direction u = surf.diffuse_reflect(p.r(), p.u(), p.current_seed());
  u /= u.norm();

  // Handle the effects of the surface albedo on the particle's weight.
  BoundaryCondition::handle_albedo(p, surf);

  p.cross_reflective_bc(surf, u);
}

//==============================================================================
// TranslationalPeriodicBC implementation
//==============================================================================

TranslationalPeriodicBC::TranslationalPeriodicBC(int i_surf, int j_surf)
  : PeriodicBC(i_surf, j_surf)
{
  Surface& surf1 {*model::surfaces[i_surf_]};
  Surface& surf2 {*model::surfaces[j_surf_]};

  // Make sure the first surface has an appropriate type.
  if (const auto* ptr = dynamic_cast<const SurfaceXPlane*>(&surf1)) {
  } else if (const auto* ptr = dynamic_cast<const SurfaceYPlane*>(&surf1)) {
  } else if (const auto* ptr = dynamic_cast<const SurfaceZPlane*>(&surf1)) {
  } else if (const auto* ptr = dynamic_cast<const SurfacePlane*>(&surf1)) {
  } else {
    throw std::invalid_argument(fmt::format(
      "Surface {} is an invalid type for "
      "translational periodic BCs. Only planes are supported for these BCs.",
      surf1.id_));
  }

  // Make sure the second surface has an appropriate type.
  if (const auto* ptr = dynamic_cast<const SurfaceXPlane*>(&surf2)) {
  } else if (const auto* ptr = dynamic_cast<const SurfaceYPlane*>(&surf2)) {
  } else if (const auto* ptr = dynamic_cast<const SurfaceZPlane*>(&surf2)) {
  } else if (const auto* ptr = dynamic_cast<const SurfacePlane*>(&surf2)) {
  } else {
    throw std::invalid_argument(fmt::format(
      "Surface {} is an invalid type for "
      "translational periodic BCs. Only planes are supported for these BCs.",
      surf2.id_));
  }

  // Compute the distance from the first surface to the origin.  Check the
  // surface evaluate function to decide if the distance is positive, negative,
  // or zero.
  Position origin {0, 0, 0};
  Direction u = surf1.normal(origin);
  double d1;
  double e1 = surf1.evaluate(origin);
  if (e1 > FP_COINCIDENT) {
    d1 = -surf1.distance(origin, -u, false);
  } else if (e1 < -FP_COINCIDENT) {
    d1 = surf1.distance(origin, u, false);
  } else {
    d1 = 0.0;
  }

  // Compute the distance from the second surface to the origin.
  double d2;
  double e2 = surf2.evaluate(origin);
  if (e2 > FP_COINCIDENT) {
    d2 = -surf2.distance(origin, -u, false);
  } else if (e2 < -FP_COINCIDENT) {
    d2 = surf2.distance(origin, u, false);
  } else {
    d2 = 0.0;
  }

  // Set the translation vector; it's length is the difference in the two
  // distances.
  translation_ = u * (d2 - d1);
}

void TranslationalPeriodicBC::handle_particle(
  Particle& p, const Surface& surf) const
{
  int i_particle_surf = p.surface_index();

  // Figure out which of the two BC surfaces were struck then find the
  // particle's new location and surface.
  Position new_r;
  int new_surface;
  if (i_particle_surf == i_surf_) {
    new_r = p.r() + translation_;
    new_surface = p.surface() > 0 ? j_surf_ + 1 : -(j_surf_ + 1);
  } else if (i_particle_surf == j_surf_) {
    new_r = p.r() - translation_;
    new_surface = p.surface() > 0 ? i_surf_ + 1 : -(i_surf_ + 1);
  } else {
    throw std::runtime_error(
      "Called BoundaryCondition::handle_particle after "
      "hitting a surface, but that surface is not recognized by the BC.");
  }

  // Handle the effects of the surface albedo on the particle's weight.
  BoundaryCondition::handle_albedo(p, surf);

  // Pass the new location and surface to the particle.
  p.cross_periodic_bc(surf, new_r, p.u(), new_surface);
}

//==============================================================================
// RotationalPeriodicBC implementation
//==============================================================================

RotationalPeriodicBC::RotationalPeriodicBC(
  int i_surf, int j_surf, PeriodicAxis axis)
  : PeriodicBC(i_surf, j_surf)
{
  Surface& surf1 {*model::surfaces[i_surf_]};
  Surface& surf2 {*model::surfaces[j_surf_]};

  // below convention for right handed coordinate system
  switch (axis) {
  case x:
    zero_axis_idx_ = 0; // x component of plane must be zero
    axis_1_idx_ = 1;    // y component independent
    axis_2_idx_ = 2;    // z component dependent
    break;
  case y:
    // for a right handed coordinate system, z should be the independent axis
    // but this would cause the y-rotation case to be different than the other
    // two. using a left handed coordinate system and a negative rotation the
    // compute angle and rotation matrix behavior mimics that of the x and z
    // cases
    zero_axis_idx_ = 1; // y component of plane must be zero
    axis_1_idx_ = 0;    // x component independent
    axis_2_idx_ = 2;    // z component dependent
    break;
  case z:
    zero_axis_idx_ = 2; // z component of plane must be zero
    axis_1_idx_ = 0;    // x component independent
    axis_2_idx_ = 1;    // y component dependent
    break;
  default:
    throw std::invalid_argument(
      fmt::format("You've specified an axis that is not x, y, or z."));
  }

  // Compute the surface normal vectors and make sure they are perpendicular
  // to the correct axis
  Direction norm1 = surf1.normal({0, 0, 0});
  Direction norm2 = surf2.normal({0, 0, 0});
  // Make sure both surfaces intersect the origin
  if (std::abs(surf1.evaluate({0, 0, 0})) > FP_COINCIDENT) {
    throw std::invalid_argument(fmt::format(
      "Rotational periodic BCs are only "
      "supported for rotations about the origin, but surface {} does not "
      "intersect the origin.",
      surf1.id_));
  }
  if (std::abs(surf2.evaluate({0, 0, 0})) > FP_COINCIDENT) {
    throw std::invalid_argument(fmt::format(
      "Rotational periodic BCs are only "
      "supported for rotations about the origin, but surface {} does not "
      "intersect the origin.",
      surf2.id_));
  }

  // reverse the angle computed if there is y-periodicity to account for
  // reversal of axes cross product direction
  angle_ = compute_periodic_rotation(norm1[axis_2_idx_], norm1[axis_1_idx_],
    norm2[axis_2_idx_], norm2[axis_1_idx_]);

  // Warn the user if the angle does not evenly divide a circle
  double rem = std::abs(std::remainder((2 * PI / angle_), 1.0));
  if (rem > FP_REL_PRECISION && rem < 1 - FP_REL_PRECISION) {
    warning(fmt::format(
      "Rotational periodic BC specified with a rotation "
      "angle of {} degrees which does not evenly divide 360 degrees.",
      angle_ * 180 / PI));
  }
}

double RotationalPeriodicBC::compute_periodic_rotation(
  double rise_1, double run_1, double rise_2, double run_2) const
{
  // Compute the BC rotation angle.  Here it is assumed that both surface
  // normal vectors point inwards---towards the valid geometry region.
  // Consequently, the rotation angle is not the difference between the two
  // normals, but is instead the difference between one normal and one
  // anti-normal.  (An incident ray on one surface must be an outgoing ray on
  // the other surface after rotation hence the anti-normal.)
  double theta1 = std::atan2(rise_1, run_1);
  double theta2 = std::atan2(rise_2, run_2) + PI;
  return theta2 - theta1;
}

void RotationalPeriodicBC::handle_particle(
  Particle& p, const Surface& surf) const
{
  int i_particle_surf = p.surface_index();

  // Figure out which of the two BC surfaces were struck to figure out if a
  // forward or backward rotation is required.  Specify the other surface as
  // the particle's new surface.
  double theta;
  int new_surface;
  if (i_particle_surf == i_surf_) {
    theta = angle_;
    new_surface = p.surface() > 0 ? -(j_surf_ + 1) : j_surf_ + 1;
  } else if (i_particle_surf == j_surf_) {
    theta = -angle_;
    new_surface = p.surface() > 0 ? -(i_surf_ + 1) : i_surf_ + 1;
  } else {
    throw std::runtime_error(
      "Called BoundaryCondition::handle_particle after "
      "hitting a surface, but that surface is not recognized by the BC.");
  }

  // Rotate the particle's position and direction about the z-axis.
  Position r = p.r();
  Direction u = p.u();
  double cos_theta = std::cos(theta);
  double sin_theta = std::sin(theta);

  Position new_r;
  new_r[zero_axis_idx_] = r[zero_axis_idx_];
  new_r[axis_1_idx_] = cos_theta * r[axis_1_idx_] - sin_theta * r[axis_2_idx_];
  new_r[axis_2_idx_] = sin_theta * r[axis_1_idx_] + cos_theta * r[axis_2_idx_];

  Direction new_u;
  new_u[zero_axis_idx_] = u[zero_axis_idx_];
  new_u[axis_1_idx_] = cos_theta * u[axis_1_idx_] - sin_theta * u[axis_2_idx_];
  new_u[axis_2_idx_] = sin_theta * u[axis_1_idx_] + cos_theta * u[axis_2_idx_];

  // Handle the effects of the surface albedo on the particle's weight.
  BoundaryCondition::handle_albedo(p, surf);

  // Pass the new location, direction, and surface to the particle.
  p.cross_periodic_bc(surf, new_r, new_u, new_surface);
}

} // namespace openmc
