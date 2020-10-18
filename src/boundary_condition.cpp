#include "openmc/boundary_condition.h"

#include <exception>

#include <fmt/core.h>

#include "openmc/constants.h"
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

  // The following blocks will resolve the type of each surface and compute the
  // appropriate translation vector

  // Check for a pair of x-planes
  if (const auto* xplane1 = dynamic_cast<const SurfaceXPlane*>(&surf1)) {
    if (const auto* xplane2 = dynamic_cast<const SurfaceXPlane*>(&surf2)) {
      translation_ = {xplane2->x0_ - xplane1->x0_, 0, 0};
    } else {
      throw std::invalid_argument(fmt::format("Invalid pair of periodic "
        "surfaces ({} and {}). For a translational periodic BC, both surfaces "
        "must be of the same type (e.g. both x-plane).", surf1.id_, surf2.id_));
    }

  // Check for a pair of y-planes
  } else if (const auto* yplane1 = dynamic_cast<const SurfaceYPlane*>(&surf1)) {
    if (const auto* yplane2 = dynamic_cast<const SurfaceYPlane*>(&surf2)) {
      translation_ = {0, yplane2->y0_ - yplane1->y0_, 0};
    } else {
      throw std::invalid_argument(fmt::format("Invalid pair of periodic "
        "surfaces ({} and {}). For a translational periodic BC, both surfaces "
        "must be of the same type (e.g. both x-plane).", surf1.id_, surf2.id_));
    }

  // Check for a pair of z-planes
  } else if (const auto* zplane1 = dynamic_cast<const SurfaceZPlane*>(&surf1)) {
    if (const auto* zplane2 = dynamic_cast<const SurfaceZPlane*>(&surf2)) {
      translation_ = {0, 0, zplane2->z0_ - zplane1->z0_};
    } else {
      throw std::invalid_argument(fmt::format("Invalid pair of periodic "
        "surfaces ({} and {}). For a translational periodic BC, both surfaces "
        "must be of the same type (e.g. both x-plane).", surf1.id_, surf2.id_));
    }
  }

  // TODO: Check for a pair of general planes

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

  // The following blocks will resolve the type of each surface and compute the
  // appropriate rotation angle

  if (const auto* xplane = dynamic_cast<const SurfaceXPlane*>(&surf1)) {
    if (const auto* yplane = dynamic_cast<const SurfaceYPlane*>(&surf2)) {
      angle_ = 0.5 * PI;
    }
  } else if (const auto* yplane = dynamic_cast<const SurfaceYPlane*>(&surf1)) {
    if (const auto* xplane = dynamic_cast<const SurfaceXPlane*>(&surf2)) {
      angle_ = -0.5 * PI;
    }
  }
}

void
RotationalPeriodicBC::handle_particle(Particle& p, const Surface& surf) const
{
  // TODO: off-by-one on surface indices throughout this function.
  int i_particle_surf = std::abs(p.surface_) - 1;

  // Figure out which of the two BC surfaces were struck to figure out if a
  // forward or backward rotation is required. Also specify the new surface.
  double theta;
  int new_surface;
  if (i_particle_surf == i_surf_) {
    theta = -angle_;
    new_surface = p.surface_ > 0 ? -(j_surf_ + 1) : j_surf_ + 1;
  } else if (i_particle_surf == j_surf_) {
    theta = angle_;
    new_surface = p.surface_ > 0 ? -(i_surf_ + 1) : i_surf_ + 1;
  } else {
    throw std::runtime_error("Called BoundaryCondition::handle_particle after "
      "hitting a surface, but that surface is not recognized by the BC.");
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
