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

RotationalPeriodicBC::RotationalPeriodicBC(int i_surf, int j_surf)
  : PeriodicBC(i_surf, j_surf)
{
  Surface& surf1 {*model::surfaces[i_surf_]};
  Surface& surf2 {*model::surfaces[j_surf_]};

  // Check the type of the first surface
  bool surf1_is_xyplane;
  if (const auto* ptr = dynamic_cast<const SurfaceXPlane*>(&surf1)) {
    surf1_is_xyplane = true;
  } else if (const auto* ptr = dynamic_cast<const SurfaceYPlane*>(&surf1)) {
    surf1_is_xyplane = true;
  } else if (const auto* ptr = dynamic_cast<const SurfacePlane*>(&surf1)) {
    surf1_is_xyplane = false;
  } else {
    throw std::invalid_argument(fmt::format(
      "Surface {} is an invalid type for "
      "rotational periodic BCs. Only x-planes, y-planes, or general planes "
      "(that are perpendicular to z) are supported for these BCs.",
      surf1.id_));
  }

  // Check the type of the second surface
  bool surf2_is_xyplane;
  if (const auto* ptr = dynamic_cast<const SurfaceXPlane*>(&surf2)) {
    surf2_is_xyplane = true;
  } else if (const auto* ptr = dynamic_cast<const SurfaceYPlane*>(&surf2)) {
    surf2_is_xyplane = true;
  } else if (const auto* ptr = dynamic_cast<const SurfacePlane*>(&surf2)) {
    surf2_is_xyplane = false;
  } else {
    throw std::invalid_argument(fmt::format(
      "Surface {} is an invalid type for "
      "rotational periodic BCs. Only x-planes, y-planes, or general planes "
      "(that are perpendicular to z) are supported for these BCs.",
      surf2.id_));
  }

  // Compute the surface normal vectors and make sure they are perpendicular
  // to the z-axis
  Direction norm1 = surf1.normal({0, 0, 0});
  Direction norm2 = surf2.normal({0, 0, 0});
  if (std::abs(norm1.z) > FP_PRECISION) {
    throw std::invalid_argument(fmt::format(
      "Rotational periodic BCs are only "
      "supported for rotations about the z-axis, but surface {} is not "
      "perpendicular to the z-axis.",
      surf1.id_));
  }
  if (std::abs(norm2.z) > FP_PRECISION) {
    throw std::invalid_argument(fmt::format(
      "Rotational periodic BCs are only "
      "supported for rotations about the z-axis, but surface {} is not "
      "perpendicular to the z-axis.",
      surf2.id_));
  }

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

  // Compute the BC rotation angle.  Here it is assumed that both surface
  // normal vectors point inwards---towards the valid geometry region.
  // Consequently, the rotation angle is not the difference between the two
  // normals, but is instead the difference between one normal and one
  // anti-normal.  (An incident ray on one surface must be an outgoing ray on
  // the other surface after rotation hence the anti-normal.)
  double theta1 = std::atan2(norm1.y, norm1.x);
  double theta2 = std::atan2(norm2.y, norm2.x) + PI;
  angle_ = theta2 - theta1;

  // Warn the user if the angle does not evenly divide a circle
  double rem = std::abs(std::remainder((2 * PI / angle_), 1.0));
  if (rem > FP_REL_PRECISION && rem < 1 - FP_REL_PRECISION) {
    warning(fmt::format(
      "Rotational periodic BC specified with a rotation "
      "angle of {} degrees which does not evenly divide 360 degrees.",
      angle_ * 180 / PI));
  }
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
  Position new_r = {
    cos_theta * r.x - sin_theta * r.y, sin_theta * r.x + cos_theta * r.y, r.z};
  Direction new_u = {
    cos_theta * u.x - sin_theta * u.y, sin_theta * u.x + cos_theta * u.y, u.z};

  // Handle the effects of the surface albedo on the particle's weight.
  BoundaryCondition::handle_albedo(p, surf);

  // Pass the new location, direction, and surface to the particle.
  p.cross_periodic_bc(surf, new_r, new_u, new_surface);
}

} // namespace openmc
