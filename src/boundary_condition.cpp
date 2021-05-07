#include "openmc/boundary_condition.h"

#include <exception>

#include <fmt/core.h>

#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/surface.h"
  
/*
  switch(type_){
    case BCType::Vacuum:                return VacuumBCzz; break;
    case BCType::Reflective:            return ReflectiveBCzz; break;
    case BCType::White:                 return WhiteBCzz; break;
    case BCType::TranslationalPeriodic: return TranslationalPeriodicBCzz; break;
    case BCType::RotationalPeriodic:    return RotationalPeriodicBCzz; break;
  }
  */

namespace openmc {

std::string BoundaryCondition::type() const
{
  switch(type_){
    case BCType::Transmission:     return "transmission";     break;
    case BCType::Vacuum:     return "vacuum";     break;
    case BCType::Reflective: return "reflective"; break;
    case BCType::White:      return "white";      break;
    default:                 return "periodic";   break;
  }
}
  
void BoundaryCondition::handle_particle(Particle& p, const Surface& surf) const
{
  switch(type_){
    case BCType::Vacuum:                return VacuumBC_handle_particle(p, surf); break;
    case BCType::Reflective:            return ReflectiveBC_handle_particle(p, surf); break;
    case BCType::White:                 return WhiteBC_handle_particle(p, surf); break;
    case BCType::TranslationalPeriodic: return TranslationalPeriodicBC_handle_particle(p, surf); break;
    case BCType::RotationalPeriodic:    return RotationalPeriodicBC_handle_particle(p, surf); break;
    default: printf("error!\n");
  }
}

BoundaryCondition::BoundaryCondition(BCType type)
{
  type_ = type;
}

BoundaryCondition::BoundaryCondition(BCType type, int i_surf, int j_surf)
{
  type_ = type;
  i_surf_ = i_surf;
  j_surf_ = j_surf;
  switch(type_)
  {
    case BCType::TranslationalPeriodic: init_TranslationalPeriodicBC(); break;
    case BCType::RotationalPeriodic:    init_RotationalPeriodicBC(); break;
    default: printf("ERROR - wrong BC type specified!\n"); break;
  }
}



//==============================================================================
// VacuumBC implementation
//==============================================================================


void
BoundaryCondition::VacuumBC_handle_particle(Particle& p, const Surface& surf) const
{
  // TODO: This does not work for some reason.
  //#pragma omp target update to(p, surf)
  //#pragma omp target
  {
    p.cross_vacuum_bc(surf);
  }
 // #pragma omp target update from(p, surf)
}


//==============================================================================
// ReflectiveBC implementation
//==============================================================================

void
BoundaryCondition::ReflectiveBC_handle_particle(Particle& p, const Surface& surf) const
{
  Direction u = surf.reflect(p.r(), p.u(), &p);
  u /= u.norm();

  p.cross_reflective_bc(surf, u);
}

//==============================================================================
// WhiteBC implementation
//==============================================================================

void
BoundaryCondition::WhiteBC_handle_particle(Particle& p, const Surface& surf) const
{
  Direction u = surf.diffuse_reflect(p.r(), p.u(), p.current_seed());
  u /= u.norm();

  p.cross_reflective_bc(surf, u);
}

//==============================================================================
// TranslationalPeriodicBC implementation
//==============================================================================

bool is_plane(Surface& surf)
{
  if(surf.type_ == Surface::SurfaceType::SurfacePlane)
    return true;
  if(surf.type_ == Surface::SurfaceType::SurfaceXPlane)
    return true;
  if(surf.type_ == Surface::SurfaceType::SurfaceYPlane)
    return true;
  if(surf.type_ == Surface::SurfaceType::SurfaceZPlane)
    return true;

  return false;
}


void BoundaryCondition::init_TranslationalPeriodicBC()
{
  Surface& surf1 {model::surfaces[i_surf_]};
  Surface& surf2 {model::surfaces[j_surf_]};

  // Make sure the first surface has an appropriate type.
  if(!is_plane(surf1))
  {
    /*
    throw std::invalid_argument(fmt::format("Surface {} is an invalid type for "
      "translational periodic BCs. Only planes are supported for these BCs.",
      surf1.id_));
      */
    printf("ERROR: surface is invalid type for BC\n");
  }

  // Make sure the second surface has an appropriate type.
  if(!is_plane(surf2))
  {
    /*
    throw std::invalid_argument(fmt::format("Surface {} is an invalid type for "
      "translational periodic BCs. Only planes are supported for these BCs.",
      surf2.id_));
      */
    printf("ERROR: surface is invalid type for BC\n");
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
BoundaryCondition::TranslationalPeriodicBC_handle_particle(Particle& p, const Surface& surf) const
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
    /*
    throw std::runtime_error("Called BoundaryCondition::handle_particle after "
      "hitting a surface, but that surface is not recognized by the BC.");
      */
    printf("ERROR- invalid BC surface!\n");
  }

  // Pass the new location and surface to the particle.
  p.cross_periodic_bc(surf, new_r, p.u(), new_surface);
}

//==============================================================================
// RotationalPeriodicBC implementation
//==============================================================================

bool is_plane_xy(Surface& surf)
{
  if(surf.type_ == Surface::SurfaceType::SurfaceXPlane)
    return true;
  if(surf.type_ == Surface::SurfaceType::SurfaceYPlane)
    return true;

  return false;
}

void BoundaryCondition::init_RotationalPeriodicBC()
{
  Surface& surf1 {model::surfaces[i_surf_]};
  Surface& surf2 {model::surfaces[j_surf_]};

  // Check the type of the first surface
  bool surf1_is_xyplane;
  if(is_plane_xy(surf1))
    surf1_is_xyplane = true;
  else if( surf1.type_ == Surface::SurfaceType::SurfacePlane )
    surf1_is_xyplane = false;
  else {
    /*
    throw std::invalid_argument(fmt::format("Surface {} is an invalid type for "
      "rotational periodic BCs. Only x-planes, y-planes, or general planes "
      "(that are perpendicular to z) are supported for these BCs.", surf1.id_));
      */
    printf("ERROR - invalid rotational BC!\n");
  }

  // Check the type of the second surface
  bool surf2_is_xyplane;
  if(is_plane_xy(surf2))
    surf1_is_xyplane = true;
  else if( surf2.type_ == Surface::SurfaceType::SurfacePlane )
    surf1_is_xyplane = false;
  else {
    /*
    throw std::invalid_argument(fmt::format("Surface {} is an invalid type for "
      "rotational periodic BCs. Only x-planes, y-planes, or general planes "
      "(that are perpendicular to z) are supported for these BCs.", surf2.id_));
      */
    printf("ERROR - invalid rotational BC!\n");
  }

  // Compute the surface normal vectors and make sure they are perpendicular
  // to the z-axis
  Direction norm1 = surf1.normal({0, 0, 0});
  Direction norm2 = surf2.normal({0, 0, 0});
  if (std::abs(norm1.z) > FP_PRECISION) {
    /*
    throw std::invalid_argument(fmt::format("Rotational periodic BCs are only "
      "supported for rotations about the z-axis, but surface {} is not "
      "perpendicular to the z-axis.", surf1.id_));
      */
    printf("ERROR - invalid rotational BC!\n");
  }
  if (std::abs(norm2.z) > FP_PRECISION) {
    /*
    throw std::invalid_argument(fmt::format("Rotational periodic BCs are only "
      "supported for rotations about the z-axis, but surface {} is not "
      "perpendicular to the z-axis.", surf2.id_));
      */
    printf("ERROR - invalid rotational BC!\n");
  }

  // Make sure both surfaces intersect the origin
  if (std::abs(surf1.evaluate({0, 0, 0})) > FP_COINCIDENT) {
    /*
    throw std::invalid_argument(fmt::format("Rotational periodic BCs are only "
      "supported for rotations about the origin, but surface {} does not "
      "intersect the origin.", surf1.id_));
      */
    printf("ERROR - invalid rotational BC!\n");
  }
  if (std::abs(surf2.evaluate({0, 0, 0})) > FP_COINCIDENT) {
    /*
    throw std::invalid_argument(fmt::format("Rotational periodic BCs are only "
      "supported for rotations about the origin, but surface {} does not "
      "intersect the origin.", surf2.id_));
      */
    printf("ERROR - invalid rotational BC!\n");
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
    /*
    warning(fmt::format("Rotational periodic BC specified with a rotation "
      "angle of {} degrees which does not evenly divide 360 degrees.",
      angle_ * 180 / PI));
      */
    printf("Warning, bad rotation angle on periodic BC!\n");
  }
}

void
BoundaryCondition::RotationalPeriodicBC_handle_particle(Particle& p, const Surface& surf) const
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
    new_surface = p.surface_ > 0 ? -(j_surf_ + 1) : j_surf_ + 1;
  } else if (i_particle_surf == j_surf_) {
    theta = -angle_;
    new_surface = p.surface_ > 0 ? -(i_surf_ + 1) : i_surf_ + 1;
  } else {
    /*
    throw std::runtime_error("Called BoundaryCondition::handle_particle after "
      "hitting a surface, but that surface is not recognized by the BC.");
      */
    printf("Error in BC handle particle!\n");
  }

  // Rotate the particle's position and direction about the z-axis.
  Position r = p.r();
  Direction u = p.u();
  double cos_theta = std::cos(theta);
  double sin_theta = std::sin(theta);
  Position new_r = {
    cos_theta*r.x - sin_theta*r.y,
    sin_theta*r.x + cos_theta*r.y,
    r.z};
  Direction new_u = {
    cos_theta*u.x - sin_theta*u.y,
    sin_theta*u.x + cos_theta*u.y,
    u.z};

  // Pass the new location, direction, and surface to the particle.
  p.cross_periodic_bc(surf, new_r, new_u, new_surface);
}

} // namespace openmc
