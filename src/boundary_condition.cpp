#include "openmc/boundary_condition.h"

#include <exception>

#include <fmt/core.h>

#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/surface.h"

namespace openmc {

//==============================================================================
// VacuumBC implementation
//==============================================================================

void
VacuumBC::handle_particle(Particle& p, const Surface& surf) const
{
  p.cross_vacuum_bc(surf);
}

//==============================================================================
// ReflectiveBC implementation
//==============================================================================

void
ReflectiveBC::handle_particle(Particle& p, const Surface& surf) const
{
  Direction u = surf.reflect(p.r(), p.u(), &p);
  u /= u.norm();

  p.cross_reflective_bc(surf, u);
}

//==============================================================================
// WhiteBC implementation
//==============================================================================

void
WhiteBC::handle_particle(Particle& p, const Surface& surf) const
{
  Direction u = surf.diffuse_reflect(p.r(), p.u(), p.current_seed());
  u /= u.norm();

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
    throw std::invalid_argument(fmt::format("Surface {} is an invalid type for "
      "translational periodic BCs. Only planes are supported for these BCs.",
      surf1.id_));
  }

  // Make sure the second surface has an appropriate type.
  if (const auto* ptr = dynamic_cast<const SurfaceXPlane*>(&surf2)) {
  } else if (const auto* ptr = dynamic_cast<const SurfaceYPlane*>(&surf2)) {
  } else if (const auto* ptr = dynamic_cast<const SurfaceZPlane*>(&surf2)) {
  } else if (const auto* ptr = dynamic_cast<const SurfacePlane*>(&surf2)) {
  } else {
    throw std::invalid_argument(fmt::format("Surface {} is an invalid type for "
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

void
TranslationalPeriodicBC::handle_particle(Particle& p, const Surface& surf) const
{
  // TODO: off-by-one on surface indices throughout this function.
  int i_particle_surf = std::abs(p.surface_) - 1;

  // Figure out which of the two BC surfaces were struck then find the
  // particle's new location and surface.
  Position new_r;
  int new_surface;
  if (i_particle_surf == i_surf_) {
    new_r = p.r() + translation_;
    new_surface = p.surface_ > 0 ? j_surf_ + 1 : -(j_surf_ + 1);
  } else if (i_particle_surf == j_surf_) {
    new_r = p.r() - translation_;
    new_surface = p.surface_ > 0 ? i_surf_ + 1 : -(i_surf_ + 1);
  } else {
    throw std::runtime_error("Called BoundaryCondition::handle_particle after "
      "hitting a surface, but that surface is not recognized by the BC.");
  }

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
    throw std::invalid_argument(fmt::format("Surface {} is an invalid type for "
      "rotational periodic BCs. Only x-planes, y-planes, or general planes "
      "(that are perpendicular to z) are supported for these BCs.", surf1.id_));
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
    throw std::invalid_argument(fmt::format("Surface {} is an invalid type for "
      "rotational periodic BCs. Only x-planes, y-planes, or general planes "
      "(that are perpendicular to z) are supported for these BCs.", surf2.id_));
  }

  // Compute the surface normal vectors and make sure they are perpendicular
  // to the z-axis
  Direction norm1 = surf1.normal({0, 0, 0});
  Direction norm2 = surf2.normal({0, 0, 0});
  if (std::abs(norm1.z) > FP_PRECISION) {
    throw std::invalid_argument(fmt::format("Rotational periodic BCs are only "
      "supported for rotations about the z-axis, but surface {} is not "
      "perpendicular to the z-axis.", surf1.id_));
  }
  if (std::abs(norm2.z) > FP_PRECISION) {
    throw std::invalid_argument(fmt::format("Rotational periodic BCs are only "
      "supported for rotations about the z-axis, but surface {} is not "
      "perpendicular to the z-axis.", surf2.id_));
  }

  // Make sure both surfaces intersect the origin
  if (std::abs(surf1.evaluate({0, 0, 0})) > FP_COINCIDENT) {
    throw std::invalid_argument(fmt::format("Rotational periodic BCs are only "
      "supported for rotations about the origin, but surface{} does not "
      "intersect the origin.", surf1.id_));
  }
  if (std::abs(surf2.evaluate({0, 0, 0})) > FP_COINCIDENT) {
    throw std::invalid_argument(fmt::format("Rotational periodic BCs are only "
      "supported for rotations about the origin, but surface{} does not "
      "intersect the origin.", surf2.id_));
  }

  // Compute the angle between the two surfaces; this is the BC rotation angle
  double theta1 = std::atan2(norm1.y, norm1.x);
  double theta2 = std::atan2(norm2.y, norm2.x);
  angle_ = theta2 - theta1;

  // Warn the user if the angle does not evenly divide a circle
  double rem = std::abs(std::remainder((2 * PI / angle_), 1.0));
  if (rem > FP_REL_PRECISION && rem < 1 - FP_REL_PRECISION) {
    warning(fmt::format("Rotational periodic BC specified with a rotation "
      "angle of {} degrees which does not evenly divide 360 degrees.",
      angle_ * 180 / PI));
  }

  // Guess whether or not the normal vectors of the two planes are aligned, i.e.
  // if an arc passing from one surface through the geometry to the other
  // surface will pass through the same "side" of both surfaces.  If the user
  // specified an x-plane and a y-plane then the geometry likely lies in the
  // first quadrant which means the normals are not aligned.  Otherwise, assume
  // the opposite.
  aligned_normals_ = !(surf1_is_xyplane && surf2_is_xyplane);
}

void
RotationalPeriodicBC::handle_particle(Particle& p, const Surface& surf) const
{
  // TODO: off-by-one on surface indices throughout this function.
  int i_particle_surf = std::abs(p.surface_) - 1;

  // Figure out which of the two BC surfaces were struck to figure out if a
  // forward or backward rotation is required.  Specify the other surface as
  // the particle's new surface.
  double theta;
  int new_surface;
  if (i_particle_surf == i_surf_) {
    theta = angle_;
    new_surface = p.surface_ > 0 ? j_surf_ + 1 : -(j_surf_ + 1);
  } else if (i_particle_surf == j_surf_) {
    theta = -angle_;
    new_surface = p.surface_ > 0 ? i_surf_ + 1 : -(i_surf_ + 1);
  } else {
    throw std::runtime_error("Called BoundaryCondition::handle_particle after "
      "hitting a surface, but that surface is not recognized by the BC.");
  }

  // If the normal vectors of the two surfaces are not aligned, then the logic
  // must be reversed for rotation and picking a new surface halfspace.
  if (not aligned_normals_) {
    theta = -theta;
    new_surface = -new_surface;
  }

  // Rotate the particle's position and direction about the z-axis.
  Position r = p.r();
  Direction u = p.u();
  Position new_r = {
    std::cos(theta)*r[0] - std::sin(theta)*r[1],
    std::sin(theta)*r[0] + std::cos(theta)*r[1],
    r[2]};
  Direction new_u = {
    std::cos(theta)*u[0] - std::sin(theta)*u[1],
    std::sin(theta)*u[0] + std::cos(theta)*u[1],
    u[2]};

  // Pass the new location, direction, and surface to the particle.
  p.cross_periodic_bc(surf, new_r, new_u, new_surface);
}

} // namespace openmc
