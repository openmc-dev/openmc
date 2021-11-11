#include "openmc/surface.h"

#include <cmath>
#include <set>
#include <utility>

#include <fmt/core.h>
#include <gsl/gsl-lite.hpp>

#include "openmc/array.h"
#include "openmc/container_util.h"
#include "openmc/error.h"
#include "openmc/hdf5_interface.h"
#include "openmc/math_functions.h"
#include "openmc/random_lcg.h"
#include "openmc/settings.h"
#include "openmc/string_utils.h"
#include "openmc/xml_interface.h"

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace model {
std::unordered_map<int, int> surface_map;
vector<unique_ptr<Surface>> surfaces;
} // namespace model

//==============================================================================
// Helper functions for reading the "coeffs" node of an XML surface element
//==============================================================================

void read_coeffs(pugi::xml_node surf_node, int surf_id, double& c1)
{
  // Check the given number of coefficients.
  std::string coeffs = get_node_value(surf_node, "coeffs");
  int n_words = word_count(coeffs);
  if (n_words != 1) {
    fatal_error(fmt::format(
      "Surface {} expects 1 coeff but was given {}", surf_id, n_words));
  }

  // Parse the coefficients.
  int stat = sscanf(coeffs.c_str(), "%lf", &c1);
  if (stat != 1) {
    fatal_error("Something went wrong reading surface coeffs");
  }
}

void read_coeffs(
  pugi::xml_node surf_node, int surf_id, double& c1, double& c2, double& c3)
{
  // Check the given number of coefficients.
  std::string coeffs = get_node_value(surf_node, "coeffs");
  int n_words = word_count(coeffs);
  if (n_words != 3) {
    fatal_error(fmt::format(
      "Surface {} expects 3 coeffs but was given {}", surf_id, n_words));
  }

  // Parse the coefficients.
  int stat = sscanf(coeffs.c_str(), "%lf %lf %lf", &c1, &c2, &c3);
  if (stat != 3) {
    fatal_error("Something went wrong reading surface coeffs");
  }
}

void read_coeffs(pugi::xml_node surf_node, int surf_id, double& c1, double& c2,
  double& c3, double& c4)
{
  // Check the given number of coefficients.
  std::string coeffs = get_node_value(surf_node, "coeffs");
  int n_words = word_count(coeffs);
  if (n_words != 4) {
    fatal_error(fmt::format(
      "Surface {} expects 4 coeffs but was given ", surf_id, n_words));
  }

  // Parse the coefficients.
  int stat = sscanf(coeffs.c_str(), "%lf %lf %lf %lf", &c1, &c2, &c3, &c4);
  if (stat != 4) {
    fatal_error("Something went wrong reading surface coeffs");
  }
}

void read_coeffs(pugi::xml_node surf_node, int surf_id, double& c1, double& c2,
  double& c3, double& c4, double& c5, double& c6, double& c7, double& c8,
  double& c9, double& c10)
{
  // Check the given number of coefficients.
  std::string coeffs = get_node_value(surf_node, "coeffs");
  int n_words = word_count(coeffs);
  if (n_words != 10) {
    fatal_error(fmt::format(
      "Surface {} expects 10 coeffs but was given {}", surf_id, n_words));
  }

  // Parse the coefficients.
  int stat = sscanf(coeffs.c_str(), "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
    &c1, &c2, &c3, &c4, &c5, &c6, &c7, &c8, &c9, &c10);
  if (stat != 10) {
    fatal_error("Something went wrong reading surface coeffs");
  }
}

//==============================================================================
// Surface implementation
//==============================================================================

Surface::Surface() {} // empty constructor

Surface::Surface(pugi::xml_node surf_node)
{
  if (check_for_node(surf_node, "id")) {
    id_ = std::stoi(get_node_value(surf_node, "id"));
    if (contains(settings::source_write_surf_id, id_)) {
      surf_source_ = true;
    }
  } else {
    fatal_error("Must specify id of surface in geometry XML file.");
  }

  if (check_for_node(surf_node, "name")) {
    name_ = get_node_value(surf_node, "name", false);
  }

  if (check_for_node(surf_node, "boundary")) {
    std::string surf_bc = get_node_value(surf_node, "boundary", true, true);

    if (surf_bc == "transmission" || surf_bc == "transmit" || surf_bc.empty()) {
      // Leave the bc_ a nullptr
    } else if (surf_bc == "vacuum") {
      bc_ = std::make_shared<VacuumBC>();
    } else if (surf_bc == "reflective" || surf_bc == "reflect" ||
               surf_bc == "reflecting") {
      bc_ = std::make_shared<ReflectiveBC>();
    } else if (surf_bc == "white") {
      bc_ = std::make_shared<WhiteBC>();
    } else if (surf_bc == "periodic") {
      // periodic BC's are handled separately
    } else {
      fatal_error(fmt::format("Unknown boundary condition \"{}\" specified "
                              "on surface {}",
        surf_bc, id_));
    }
  }
}

bool Surface::sense(Position r, Direction u) const
{
  // Evaluate the surface equation at the particle's coordinates to determine
  // which side the particle is on.
  const double f = evaluate(r);

  // Check which side of surface the point is on.
  if (std::abs(f) < FP_COINCIDENT) {
    // Particle may be coincident with this surface. To determine the sense, we
    // look at the direction of the particle relative to the surface normal (by
    // default in the positive direction) via their dot product.
    return u.dot(normal(r)) > 0.0;
  }
  return f > 0.0;
}

Direction Surface::reflect(Position r, Direction u, Particle* p) const
{
  // Determine projection of direction onto normal and squared magnitude of
  // normal.
  Direction n = normal(r);

  // Reflect direction according to normal.
  return u.reflect(n);
}

Direction Surface::diffuse_reflect(
  Position r, Direction u, uint64_t* seed) const
{
  // Diffuse reflect direction according to the normal.
  // cosine distribution

  Direction n = this->normal(r);
  n /= n.norm();
  const double projection = n.dot(u);

  // sample from inverse function, u=sqrt(rand) since p(u)=2u, so F(u)=u^2
  const double mu =
    (projection >= 0.0) ? -std::sqrt(prn(seed)) : std::sqrt(prn(seed));

  // sample azimuthal distribution uniformly
  u = rotate_angle(n, mu, nullptr, seed);

  // normalize the direction
  return u / u.norm();
}

void Surface::to_hdf5(hid_t group_id) const
{
  hid_t surf_group = create_group(group_id, fmt::format("surface {}", id_));

  if (geom_type_ == GeometryType::DAG) {
    write_string(surf_group, "geom_type", "dagmc", false);
  } else if (geom_type_ == GeometryType::CSG) {
    write_string(surf_group, "geom_type", "csg", false);

    if (bc_) {
      write_string(surf_group, "boundary_type", bc_->type(), false);
    } else {
      write_string(surf_group, "boundary_type", "transmission", false);
    }
  }

  if (!name_.empty()) {
    write_string(surf_group, "name", name_, false);
  }

  to_hdf5_inner(surf_group);

  close_group(surf_group);
}

CSGSurface::CSGSurface() : Surface {}
{
  geom_type_ = GeometryType::CSG;
};
CSGSurface::CSGSurface(pugi::xml_node surf_node) : Surface {surf_node}
{
  geom_type_ = GeometryType::CSG;
};

//==============================================================================
// Generic functions for x-, y-, and z-, planes.
//==============================================================================

// The template parameter indicates the axis normal to the plane.
template<int i>
double axis_aligned_plane_distance(
  Position r, Direction u, bool coincident, double offset)
{
  const double f = offset - r[i];
  if (coincident || std::abs(f) < FP_COINCIDENT || u[i] == 0.0)
    return INFTY;
  const double d = f / u[i];
  if (d < 0.0)
    return INFTY;
  return d;
}

//==============================================================================
// SurfaceXPlane implementation
//==============================================================================

SurfaceXPlane::SurfaceXPlane(pugi::xml_node surf_node) : CSGSurface(surf_node)
{
  read_coeffs(surf_node, id_, x0_);
}

double SurfaceXPlane::evaluate(Position r) const
{
  return r.x - x0_;
}

double SurfaceXPlane::distance(Position r, Direction u, bool coincident) const
{
  return axis_aligned_plane_distance<0>(r, u, coincident, x0_);
}

Direction SurfaceXPlane::normal(Position r) const
{
  return {1., 0., 0.};
}

void SurfaceXPlane::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "x-plane", false);
  array<double, 1> coeffs {{x0_}};
  write_dataset(group_id, "coefficients", coeffs);
}

BoundingBox SurfaceXPlane::bounding_box(bool pos_side) const
{
  if (pos_side) {
    return {x0_, INFTY, -INFTY, INFTY, -INFTY, INFTY};
  } else {
    return {-INFTY, x0_, -INFTY, INFTY, -INFTY, INFTY};
  }
}

//==============================================================================
// SurfaceYPlane implementation
//==============================================================================

SurfaceYPlane::SurfaceYPlane(pugi::xml_node surf_node) : CSGSurface(surf_node)
{
  read_coeffs(surf_node, id_, y0_);
}

double SurfaceYPlane::evaluate(Position r) const
{
  return r.y - y0_;
}

double SurfaceYPlane::distance(Position r, Direction u, bool coincident) const
{
  return axis_aligned_plane_distance<1>(r, u, coincident, y0_);
}

Direction SurfaceYPlane::normal(Position r) const
{
  return {0., 1., 0.};
}

void SurfaceYPlane::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "y-plane", false);
  array<double, 1> coeffs {{y0_}};
  write_dataset(group_id, "coefficients", coeffs);
}

BoundingBox SurfaceYPlane::bounding_box(bool pos_side) const
{
  if (pos_side) {
    return {-INFTY, INFTY, y0_, INFTY, -INFTY, INFTY};
  } else {
    return {-INFTY, INFTY, -INFTY, y0_, -INFTY, INFTY};
  }
}

//==============================================================================
// SurfaceZPlane implementation
//==============================================================================

SurfaceZPlane::SurfaceZPlane(pugi::xml_node surf_node) : CSGSurface(surf_node)
{
  read_coeffs(surf_node, id_, z0_);
}

double SurfaceZPlane::evaluate(Position r) const
{
  return r.z - z0_;
}

double SurfaceZPlane::distance(Position r, Direction u, bool coincident) const
{
  return axis_aligned_plane_distance<2>(r, u, coincident, z0_);
}

Direction SurfaceZPlane::normal(Position r) const
{
  return {0., 0., 1.};
}

void SurfaceZPlane::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "z-plane", false);
  array<double, 1> coeffs {{z0_}};
  write_dataset(group_id, "coefficients", coeffs);
}

BoundingBox SurfaceZPlane::bounding_box(bool pos_side) const
{
  if (pos_side) {
    return {-INFTY, INFTY, -INFTY, INFTY, z0_, INFTY};
  } else {
    return {-INFTY, INFTY, -INFTY, INFTY, -INFTY, z0_};
  }
}

//==============================================================================
// SurfacePlane implementation
//==============================================================================

SurfacePlane::SurfacePlane(pugi::xml_node surf_node) : CSGSurface(surf_node)
{
  read_coeffs(surf_node, id_, A_, B_, C_, D_);
}

double SurfacePlane::evaluate(Position r) const
{
  return A_ * r.x + B_ * r.y + C_ * r.z - D_;
}

double SurfacePlane::distance(Position r, Direction u, bool coincident) const
{
  const double f = A_ * r.x + B_ * r.y + C_ * r.z - D_;
  const double projection = A_ * u.x + B_ * u.y + C_ * u.z;
  if (coincident || std::abs(f) < FP_COINCIDENT || projection == 0.0) {
    return INFTY;
  } else {
    const double d = -f / projection;
    if (d < 0.0)
      return INFTY;
    return d;
  }
}

Direction SurfacePlane::normal(Position r) const
{
  return {A_, B_, C_};
}

void SurfacePlane::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "plane", false);
  array<double, 4> coeffs {{A_, B_, C_, D_}};
  write_dataset(group_id, "coefficients", coeffs);
}

//==============================================================================
// Generic functions for x-, y-, and z-, cylinders
//==============================================================================

// The template parameters indicate the axes perpendicular to the axis of the
// cylinder.  offset1 and offset2 should correspond with i1 and i2,
// respectively.
template<int i1, int i2>
double axis_aligned_cylinder_evaluate(
  Position r, double offset1, double offset2, double radius)
{
  const double r1 = r.get<i1>() - offset1;
  const double r2 = r.get<i2>() - offset2;
  return r1 * r1 + r2 * r2 - radius * radius;
}

// The first template parameter indicates which axis the cylinder is aligned to.
// The other two parameters indicate the other two axes.  offset1 and offset2
// should correspond with i2 and i3, respectively.
template<int i1, int i2, int i3>
double axis_aligned_cylinder_distance(Position r, Direction u, bool coincident,
  double offset1, double offset2, double radius)
{
  const double a = 1.0 - u.get<i1>() * u.get<i1>(); // u^2 + v^2
  if (a == 0.0)
    return INFTY;

  const double r2 = r.get<i2>() - offset1;
  const double r3 = r.get<i3>() - offset2;
  const double k = r2 * u.get<i2>() + r3 * u.get<i3>();
  const double c = r2 * r2 + r3 * r3 - radius * radius;
  const double quad = k * k - a * c;

  if (quad < 0.0) {
    // No intersection with cylinder.
    return INFTY;

  } else if (coincident || std::abs(c) < FP_COINCIDENT) {
    // Particle is on the cylinder, thus one distance is positive/negative
    // and the other is zero. The sign of k determines if we are facing in or
    // out.
    if (k >= 0.0) {
      return INFTY;
    } else {
      return (-k + sqrt(quad)) / a;
    }

  } else if (c < 0.0) {
    // Particle is inside the cylinder, thus one distance must be negative
    // and one must be positive. The positive distance will be the one with
    // negative sign on sqrt(quad).
    return (-k + sqrt(quad)) / a;

  } else {
    // Particle is outside the cylinder, thus both distances are either
    // positive or negative. If positive, the smaller distance is the one
    // with positive sign on sqrt(quad).
    const double d = (-k - sqrt(quad)) / a;
    if (d < 0.0)
      return INFTY;
    return d;
  }
}

// The first template parameter indicates which axis the cylinder is aligned to.
// The other two parameters indicate the other two axes.  offset1 and offset2
// should correspond with i2 and i3, respectively.
template<int i1, int i2, int i3>
Direction axis_aligned_cylinder_normal(
  Position r, double offset1, double offset2)
{
  Direction u;
  u.get<i2>() = 2.0 * (r.get<i2>() - offset1);
  u.get<i3>() = 2.0 * (r.get<i3>() - offset2);
  u.get<i1>() = 0.0;
  return u;
}

//==============================================================================
// SurfaceXCylinder implementation
//==============================================================================

SurfaceXCylinder::SurfaceXCylinder(pugi::xml_node surf_node)
  : CSGSurface(surf_node)
{
  read_coeffs(surf_node, id_, y0_, z0_, radius_);
}

double SurfaceXCylinder::evaluate(Position r) const
{
  return axis_aligned_cylinder_evaluate<1, 2>(r, y0_, z0_, radius_);
}

double SurfaceXCylinder::distance(
  Position r, Direction u, bool coincident) const
{
  return axis_aligned_cylinder_distance<0, 1, 2>(
    r, u, coincident, y0_, z0_, radius_);
}

Direction SurfaceXCylinder::normal(Position r) const
{
  return axis_aligned_cylinder_normal<0, 1, 2>(r, y0_, z0_);
}

void SurfaceXCylinder::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "x-cylinder", false);
  array<double, 3> coeffs {{y0_, z0_, radius_}};
  write_dataset(group_id, "coefficients", coeffs);
}

BoundingBox SurfaceXCylinder::bounding_box(bool pos_side) const
{
  if (!pos_side) {
    return {-INFTY, INFTY, y0_ - radius_, y0_ + radius_, z0_ - radius_,
      z0_ + radius_};
  } else {
    return {};
  }
}
//==============================================================================
// SurfaceYCylinder implementation
//==============================================================================

SurfaceYCylinder::SurfaceYCylinder(pugi::xml_node surf_node)
  : CSGSurface(surf_node)
{
  read_coeffs(surf_node, id_, x0_, z0_, radius_);
}

double SurfaceYCylinder::evaluate(Position r) const
{
  return axis_aligned_cylinder_evaluate<0, 2>(r, x0_, z0_, radius_);
}

double SurfaceYCylinder::distance(
  Position r, Direction u, bool coincident) const
{
  return axis_aligned_cylinder_distance<1, 0, 2>(
    r, u, coincident, x0_, z0_, radius_);
}

Direction SurfaceYCylinder::normal(Position r) const
{
  return axis_aligned_cylinder_normal<1, 0, 2>(r, x0_, z0_);
}

void SurfaceYCylinder::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "y-cylinder", false);
  array<double, 3> coeffs {{x0_, z0_, radius_}};
  write_dataset(group_id, "coefficients", coeffs);
}

BoundingBox SurfaceYCylinder::bounding_box(bool pos_side) const
{
  if (!pos_side) {
    return {x0_ - radius_, x0_ + radius_, -INFTY, INFTY, z0_ - radius_,
      z0_ + radius_};
  } else {
    return {};
  }
}

//==============================================================================
// SurfaceZCylinder implementation
//==============================================================================

SurfaceZCylinder::SurfaceZCylinder(pugi::xml_node surf_node)
  : CSGSurface(surf_node)
{
  read_coeffs(surf_node, id_, x0_, y0_, radius_);
}

double SurfaceZCylinder::evaluate(Position r) const
{
  return axis_aligned_cylinder_evaluate<0, 1>(r, x0_, y0_, radius_);
}

double SurfaceZCylinder::distance(
  Position r, Direction u, bool coincident) const
{
  return axis_aligned_cylinder_distance<2, 0, 1>(
    r, u, coincident, x0_, y0_, radius_);
}

Direction SurfaceZCylinder::normal(Position r) const
{
  return axis_aligned_cylinder_normal<2, 0, 1>(r, x0_, y0_);
}

void SurfaceZCylinder::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "z-cylinder", false);
  array<double, 3> coeffs {{x0_, y0_, radius_}};
  write_dataset(group_id, "coefficients", coeffs);
}

BoundingBox SurfaceZCylinder::bounding_box(bool pos_side) const
{
  if (!pos_side) {
    return {x0_ - radius_, x0_ + radius_, y0_ - radius_, y0_ + radius_, -INFTY,
      INFTY};
  } else {
    return {};
  }
}

//==============================================================================
// SurfaceSphere implementation
//==============================================================================

SurfaceSphere::SurfaceSphere(pugi::xml_node surf_node) : CSGSurface(surf_node)
{
  read_coeffs(surf_node, id_, x0_, y0_, z0_, radius_);
}

double SurfaceSphere::evaluate(Position r) const
{
  const double x = r.x - x0_;
  const double y = r.y - y0_;
  const double z = r.z - z0_;
  return x * x + y * y + z * z - radius_ * radius_;
}

double SurfaceSphere::distance(Position r, Direction u, bool coincident) const
{
  const double x = r.x - x0_;
  const double y = r.y - y0_;
  const double z = r.z - z0_;
  const double k = x * u.x + y * u.y + z * u.z;
  const double c = x * x + y * y + z * z - radius_ * radius_;
  const double quad = k * k - c;

  if (quad < 0.0) {
    // No intersection with sphere.
    return INFTY;

  } else if (coincident || std::abs(c) < FP_COINCIDENT) {
    // Particle is on the sphere, thus one distance is positive/negative and
    // the other is zero. The sign of k determines if we are facing in or out.
    if (k >= 0.0) {
      return INFTY;
    } else {
      return -k + sqrt(quad);
    }

  } else if (c < 0.0) {
    // Particle is inside the sphere, thus one distance must be negative and
    // one must be positive. The positive distance will be the one with
    // negative sign on sqrt(quad)
    return -k + sqrt(quad);

  } else {
    // Particle is outside the sphere, thus both distances are either positive
    // or negative. If positive, the smaller distance is the one with positive
    // sign on sqrt(quad).
    const double d = -k - sqrt(quad);
    if (d < 0.0)
      return INFTY;
    return d;
  }
}

Direction SurfaceSphere::normal(Position r) const
{
  return {2.0 * (r.x - x0_), 2.0 * (r.y - y0_), 2.0 * (r.z - z0_)};
}

void SurfaceSphere::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "sphere", false);
  array<double, 4> coeffs {{x0_, y0_, z0_, radius_}};
  write_dataset(group_id, "coefficients", coeffs);
}

BoundingBox SurfaceSphere::bounding_box(bool pos_side) const
{
  if (!pos_side) {
    return {x0_ - radius_, x0_ + radius_, y0_ - radius_, y0_ + radius_,
      z0_ - radius_, z0_ + radius_};
  } else {
    return {};
  }
}

//==============================================================================
// Generic functions for x-, y-, and z-, cones
//==============================================================================

// The first template parameter indicates which axis the cone is aligned to.
// The other two parameters indicate the other two axes.  offset1, offset2,
// and offset3 should correspond with i1, i2, and i3, respectively.
template<int i1, int i2, int i3>
double axis_aligned_cone_evaluate(
  Position r, double offset1, double offset2, double offset3, double radius_sq)
{
  const double r1 = r.get<i1>() - offset1;
  const double r2 = r.get<i2>() - offset2;
  const double r3 = r.get<i3>() - offset3;
  return r2 * r2 + r3 * r3 - radius_sq * r1 * r1;
}

// The first template parameter indicates which axis the cone is aligned to.
// The other two parameters indicate the other two axes.  offset1, offset2,
// and offset3 should correspond with i1, i2, and i3, respectively.
template<int i1, int i2, int i3>
double axis_aligned_cone_distance(Position r, Direction u, bool coincident,
  double offset1, double offset2, double offset3, double radius_sq)
{
  const double r1 = r.get<i1>() - offset1;
  const double r2 = r.get<i2>() - offset2;
  const double r3 = r.get<i3>() - offset3;
  const double a = u.get<i2>() * u.get<i2>() + u.get<i3>() * u.get<i3>() -
                   radius_sq * u.get<i1>() * u.get<i1>();
  const double k =
    r2 * u.get<i2>() + r3 * u.get<i3>() - radius_sq * r1 * u.get<i1>();
  const double c = r2 * r2 + r3 * r3 - radius_sq * r1 * r1;
  double quad = k * k - a * c;

  double d;

  if (quad < 0.0) {
    // No intersection with cone.
    return INFTY;

  } else if (coincident || std::abs(c) < FP_COINCIDENT) {
    // Particle is on the cone, thus one distance is positive/negative
    // and the other is zero. The sign of k determines if we are facing in or
    // out.
    if (k >= 0.0) {
      d = (-k - sqrt(quad)) / a;
    } else {
      d = (-k + sqrt(quad)) / a;
    }

  } else {
    // Calculate both solutions to the quadratic.
    quad = sqrt(quad);
    d = (-k - quad) / a;
    const double b = (-k + quad) / a;

    // Determine the smallest positive solution.
    if (d < 0.0) {
      if (b > 0.0)
        d = b;
    } else {
      if (b > 0.0) {
        if (b < d)
          d = b;
      }
    }
  }

  // If the distance was negative, set boundary distance to infinity.
  if (d <= 0.0)
    return INFTY;
  return d;
}

// The first template parameter indicates which axis the cone is aligned to.
// The other two parameters indicate the other two axes.  offset1, offset2,
// and offset3 should correspond with i1, i2, and i3, respectively.
template<int i1, int i2, int i3>
Direction axis_aligned_cone_normal(
  Position r, double offset1, double offset2, double offset3, double radius_sq)
{
  Direction u;
  u.get<i1>() = -2.0 * radius_sq * (r.get<i1>() - offset1);
  u.get<i2>() = 2.0 * (r.get<i2>() - offset2);
  u.get<i3>() = 2.0 * (r.get<i3>() - offset3);
  return u;
}

//==============================================================================
// SurfaceXCone implementation
//==============================================================================

SurfaceXCone::SurfaceXCone(pugi::xml_node surf_node) : CSGSurface(surf_node)
{
  read_coeffs(surf_node, id_, x0_, y0_, z0_, radius_sq_);
}

double SurfaceXCone::evaluate(Position r) const
{
  return axis_aligned_cone_evaluate<0, 1, 2>(r, x0_, y0_, z0_, radius_sq_);
}

double SurfaceXCone::distance(Position r, Direction u, bool coincident) const
{
  return axis_aligned_cone_distance<0, 1, 2>(
    r, u, coincident, x0_, y0_, z0_, radius_sq_);
}

Direction SurfaceXCone::normal(Position r) const
{
  return axis_aligned_cone_normal<0, 1, 2>(r, x0_, y0_, z0_, radius_sq_);
}

void SurfaceXCone::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "x-cone", false);
  array<double, 4> coeffs {{x0_, y0_, z0_, radius_sq_}};
  write_dataset(group_id, "coefficients", coeffs);
}

//==============================================================================
// SurfaceYCone implementation
//==============================================================================

SurfaceYCone::SurfaceYCone(pugi::xml_node surf_node) : CSGSurface(surf_node)
{
  read_coeffs(surf_node, id_, x0_, y0_, z0_, radius_sq_);
}

double SurfaceYCone::evaluate(Position r) const
{
  return axis_aligned_cone_evaluate<1, 0, 2>(r, y0_, x0_, z0_, radius_sq_);
}

double SurfaceYCone::distance(Position r, Direction u, bool coincident) const
{
  return axis_aligned_cone_distance<1, 0, 2>(
    r, u, coincident, y0_, x0_, z0_, radius_sq_);
}

Direction SurfaceYCone::normal(Position r) const
{
  return axis_aligned_cone_normal<1, 0, 2>(r, y0_, x0_, z0_, radius_sq_);
}

void SurfaceYCone::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "y-cone", false);
  array<double, 4> coeffs {{x0_, y0_, z0_, radius_sq_}};
  write_dataset(group_id, "coefficients", coeffs);
}

//==============================================================================
// SurfaceZCone implementation
//==============================================================================

SurfaceZCone::SurfaceZCone(pugi::xml_node surf_node) : CSGSurface(surf_node)
{
  read_coeffs(surf_node, id_, x0_, y0_, z0_, radius_sq_);
}

double SurfaceZCone::evaluate(Position r) const
{
  return axis_aligned_cone_evaluate<2, 0, 1>(r, z0_, x0_, y0_, radius_sq_);
}

double SurfaceZCone::distance(Position r, Direction u, bool coincident) const
{
  return axis_aligned_cone_distance<2, 0, 1>(
    r, u, coincident, z0_, x0_, y0_, radius_sq_);
}

Direction SurfaceZCone::normal(Position r) const
{
  return axis_aligned_cone_normal<2, 0, 1>(r, z0_, x0_, y0_, radius_sq_);
}

void SurfaceZCone::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "z-cone", false);
  array<double, 4> coeffs {{x0_, y0_, z0_, radius_sq_}};
  write_dataset(group_id, "coefficients", coeffs);
}

//==============================================================================
// SurfaceQuadric implementation
//==============================================================================

SurfaceQuadric::SurfaceQuadric(pugi::xml_node surf_node) : CSGSurface(surf_node)
{
  read_coeffs(surf_node, id_, A_, B_, C_, D_, E_, F_, G_, H_, J_, K_);
}

double SurfaceQuadric::evaluate(Position r) const
{
  const double x = r.x;
  const double y = r.y;
  const double z = r.z;
  return x * (A_ * x + D_ * y + G_) + y * (B_ * y + E_ * z + H_) +
         z * (C_ * z + F_ * x + J_) + K_;
}

double SurfaceQuadric::distance(
  Position r, Direction ang, bool coincident) const
{
  const double& x = r.x;
  const double& y = r.y;
  const double& z = r.z;
  const double& u = ang.x;
  const double& v = ang.y;
  const double& w = ang.z;

  const double a =
    A_ * u * u + B_ * v * v + C_ * w * w + D_ * u * v + E_ * v * w + F_ * u * w;
  const double k = A_ * u * x + B_ * v * y + C_ * w * z +
                   0.5 * (D_ * (u * y + v * x) + E_ * (v * z + w * y) +
                           F_ * (w * x + u * z) + G_ * u + H_ * v + J_ * w);
  const double c = A_ * x * x + B_ * y * y + C_ * z * z + D_ * x * y +
                   E_ * y * z + F_ * x * z + G_ * x + H_ * y + J_ * z + K_;
  double quad = k * k - a * c;

  double d;

  if (quad < 0.0) {
    // No intersection with surface.
    return INFTY;

  } else if (coincident || std::abs(c) < FP_COINCIDENT) {
    // Particle is on the surface, thus one distance is positive/negative and
    // the other is zero. The sign of k determines which distance is zero and
    // which is not. Additionally, if a is zero, it means the particle is on
    // a plane-like surface.
    if (a == 0.0) {
      d = INFTY; // see the below explanation
    } else if (k >= 0.0) {
      d = (-k - sqrt(quad)) / a;
    } else {
      d = (-k + sqrt(quad)) / a;
    }

  } else if (a == 0.0) {
    // Given the orientation of the particle, the quadric looks like a plane in
    // this case, and thus we have only one solution despite potentially having
    // quad > 0.0. While the term under the square root may be real, in one
    // case of the +/- of the quadratic formula, 0/0 results, and in another, a
    // finite value over 0 results. Applying L'Hopital's to the 0/0 case gives
    // the below. Alternatively this can be found by simply putting a=0 in the
    // equation ax^2 + bx + c = 0.
    d = -0.5 * c / k;
  } else {
    // Calculate both solutions to the quadratic.
    quad = sqrt(quad);
    d = (-k - quad) / a;
    double b = (-k + quad) / a;

    // Determine the smallest positive solution.
    if (d < 0.0) {
      if (b > 0.0)
        d = b;
    } else {
      if (b > 0.0) {
        if (b < d)
          d = b;
      }
    }
  }

  // If the distance was negative, set boundary distance to infinity.
  if (d <= 0.0)
    return INFTY;
  return d;
}

Direction SurfaceQuadric::normal(Position r) const
{
  const double& x = r.x;
  const double& y = r.y;
  const double& z = r.z;
  return {2.0 * A_ * x + D_ * y + F_ * z + G_,
    2.0 * B_ * y + D_ * x + E_ * z + H_, 2.0 * C_ * z + E_ * y + F_ * x + J_};
}

void SurfaceQuadric::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "quadric", false);
  array<double, 10> coeffs {{A_, B_, C_, D_, E_, F_, G_, H_, J_, K_}};
  write_dataset(group_id, "coefficients", coeffs);
}

//==============================================================================

void read_surfaces(pugi::xml_node node)
{
  // Count the number of surfaces
  int n_surfaces = 0;
  for (pugi::xml_node surf_node : node.children("surface")) {
    n_surfaces++;
  }

  // Loop over XML surface elements and populate the array.  Keep track of
  // periodic surfaces.
  model::surfaces.reserve(n_surfaces);
  std::set<std::pair<int, int>> periodic_pairs;
  {
    pugi::xml_node surf_node;
    int i_surf;
    for (surf_node = node.child("surface"), i_surf = 0; surf_node;
         surf_node = surf_node.next_sibling("surface"), i_surf++) {
      std::string surf_type = get_node_value(surf_node, "type", true, true);

      // Allocate and initialize the new surface

      if (surf_type == "x-plane") {
        model::surfaces.push_back(make_unique<SurfaceXPlane>(surf_node));

      } else if (surf_type == "y-plane") {
        model::surfaces.push_back(make_unique<SurfaceYPlane>(surf_node));

      } else if (surf_type == "z-plane") {
        model::surfaces.push_back(make_unique<SurfaceZPlane>(surf_node));

      } else if (surf_type == "plane") {
        model::surfaces.push_back(make_unique<SurfacePlane>(surf_node));

      } else if (surf_type == "x-cylinder") {
        model::surfaces.push_back(make_unique<SurfaceXCylinder>(surf_node));

      } else if (surf_type == "y-cylinder") {
        model::surfaces.push_back(make_unique<SurfaceYCylinder>(surf_node));

      } else if (surf_type == "z-cylinder") {
        model::surfaces.push_back(make_unique<SurfaceZCylinder>(surf_node));

      } else if (surf_type == "sphere") {
        model::surfaces.push_back(make_unique<SurfaceSphere>(surf_node));

      } else if (surf_type == "x-cone") {
        model::surfaces.push_back(make_unique<SurfaceXCone>(surf_node));

      } else if (surf_type == "y-cone") {
        model::surfaces.push_back(make_unique<SurfaceYCone>(surf_node));

      } else if (surf_type == "z-cone") {
        model::surfaces.push_back(make_unique<SurfaceZCone>(surf_node));

      } else if (surf_type == "quadric") {
        model::surfaces.push_back(make_unique<SurfaceQuadric>(surf_node));

      } else {
        fatal_error(fmt::format("Invalid surface type, \"{}\"", surf_type));
      }

      // Check for a periodic surface
      if (check_for_node(surf_node, "boundary")) {
        std::string surf_bc = get_node_value(surf_node, "boundary", true, true);
        if (surf_bc == "periodic") {
          if (check_for_node(surf_node, "periodic_surface_id")) {
            int i_periodic =
              std::stoi(get_node_value(surf_node, "periodic_surface_id"));
            int lo_id = std::min(model::surfaces.back()->id_, i_periodic);
            int hi_id = std::max(model::surfaces.back()->id_, i_periodic);
            periodic_pairs.insert({lo_id, hi_id});
          } else {
            periodic_pairs.insert({model::surfaces.back()->id_, -1});
          }
        }
      }
    }
  }

  // Fill the surface map
  for (int i_surf = 0; i_surf < model::surfaces.size(); i_surf++) {
    int id = model::surfaces[i_surf]->id_;
    auto in_map = model::surface_map.find(id);
    if (in_map == model::surface_map.end()) {
      model::surface_map[id] = i_surf;
    } else {
      fatal_error(
        fmt::format("Two or more surfaces use the same unique ID: {}", id));
    }
  }

  // Resolve unpaired periodic surfaces.  A lambda function is used with
  // std::find_if to identify the unpaired surfaces.
  auto is_unresolved_pair = [](const std::pair<int, int> p) {
    return p.second == -1;
  };
  auto first_unresolved = std::find_if(
    periodic_pairs.begin(), periodic_pairs.end(), is_unresolved_pair);
  if (first_unresolved != periodic_pairs.end()) {
    // Found one unpaired surface; search for a second one
    auto next_elem = first_unresolved;
    next_elem++;
    auto second_unresolved =
      std::find_if(next_elem, periodic_pairs.end(), is_unresolved_pair);
    if (second_unresolved == periodic_pairs.end()) {
      fatal_error("Found only one periodic surface without a specified partner."
                  " Please specify the partner for each periodic surface.");
    }

    // Make sure there isn't a third unpaired surface
    next_elem = second_unresolved;
    next_elem++;
    auto third_unresolved =
      std::find_if(next_elem, periodic_pairs.end(), is_unresolved_pair);
    if (third_unresolved != periodic_pairs.end()) {
      fatal_error(
        "Found at least three periodic surfaces without a specified "
        "partner. Please specify the partner for each periodic surface.");
    }

    // Add the completed pair and remove the old, unpaired entries
    int lo_id = std::min(first_unresolved->first, second_unresolved->first);
    int hi_id = std::max(first_unresolved->first, second_unresolved->first);
    periodic_pairs.insert({lo_id, hi_id});
    periodic_pairs.erase(first_unresolved);
    periodic_pairs.erase(second_unresolved);
  }

  // Assign the periodic boundary conditions
  for (auto periodic_pair : periodic_pairs) {
    int i_surf = model::surface_map[periodic_pair.first];
    int j_surf = model::surface_map[periodic_pair.second];
    Surface& surf1 {*model::surfaces[i_surf]};
    Surface& surf2 {*model::surfaces[j_surf]};

    // Compute the dot product of the surface normals
    Direction norm1 = surf1.normal({0, 0, 0});
    Direction norm2 = surf2.normal({0, 0, 0});
    norm1 /= norm1.norm();
    norm2 /= norm2.norm();
    double dot_prod = norm1.dot(norm2);

    // If the dot product is 1 (to within floating point precision) then the
    // planes are parallel which indicates a translational periodic boundary
    // condition.  Otherwise, it is a rotational periodic BC.
    if (std::abs(1.0 - dot_prod) < FP_PRECISION) {
      surf1.bc_ = std::make_shared<TranslationalPeriodicBC>(i_surf, j_surf);
      surf2.bc_ = surf1.bc_;
    } else {
      surf1.bc_ = std::make_shared<RotationalPeriodicBC>(i_surf, j_surf);
      surf2.bc_ = surf1.bc_;
    }
  }
}

void free_memory_surfaces()
{
  model::surfaces.clear();
  model::surface_map.clear();
}

} // namespace openmc
