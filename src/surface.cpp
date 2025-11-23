#include "openmc/surface.h"

#include <cmath>
#include <complex>
#include <initializer_list>
#include <set>
#include <utility>

#include <fmt/core.h>

#include "openmc/array.h"
#include "openmc/container_util.h"
#include "openmc/error.h"
#include "openmc/external/quartic_solver.h"
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

void read_coeffs(
  pugi::xml_node surf_node, int surf_id, std::initializer_list<double*> coeffs)
{
  // Check the given number of coefficients.
  auto coeffs_file = get_node_array<double>(surf_node, "coeffs");
  if (coeffs_file.size() != coeffs.size()) {
    fatal_error(
      fmt::format("Surface {} expects {} coefficient but was given {}", surf_id,
        coeffs.size(), coeffs_file.size()));
  }

  // Copy the coefficients
  int i = 0;
  for (auto c : coeffs) {
    *c = coeffs_file[i++];
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
    if (contains(settings::source_write_surf_id, id_) ||
        settings::source_write_surf_id.empty()) {
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
      bc_ = make_unique<VacuumBC>();
    } else if (surf_bc == "reflective" || surf_bc == "reflect" ||
               surf_bc == "reflecting") {
      bc_ = make_unique<ReflectiveBC>();
    } else if (surf_bc == "white") {
      bc_ = make_unique<WhiteBC>();
    } else if (surf_bc == "periodic") {
      // Periodic BCs are handled separately
    } else {
      fatal_error(fmt::format("Unknown boundary condition \"{}\" specified "
                              "on surface {}",
        surf_bc, id_));
    }

    if (check_for_node(surf_node, "albedo") && bc_) {
      double surf_alb = std::stod(get_node_value(surf_node, "albedo"));

      if (surf_alb < 0.0)
        fatal_error(fmt::format("Surface {} has an albedo of {}. "
                                "Albedo values must be positive.",
          id_, surf_alb));

      if (surf_alb > 1.0)
        warning(fmt::format("Surface {} has an albedo of {}. "
                            "Albedos greater than 1 may cause "
                            "unphysical behaviour.",
          id_, surf_alb));

      bc_->set_albedo(surf_alb);
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

Direction Surface::reflect(Position r, Direction u, GeometryState* p) const
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

  if (geom_type() == GeometryType::DAG) {
    write_string(surf_group, "geom_type", "dagmc", false);
  } else if (geom_type() == GeometryType::CSG) {
    write_string(surf_group, "geom_type", "csg", false);

    if (bc_) {
      write_string(surf_group, "boundary_type", bc_->type(), false);
      bc_->to_hdf5(surf_group);

      // write periodic surface ID
      if (bc_->type() == "periodic") {
        auto pbc = dynamic_cast<PeriodicBC*>(bc_.get());
        Surface& surf1 {*model::surfaces[pbc->i_surf()]};
        Surface& surf2 {*model::surfaces[pbc->j_surf()]};

        if (id_ == surf1.id_) {
          write_dataset(surf_group, "periodic_surface_id", surf2.id_);
        } else {
          write_dataset(surf_group, "periodic_surface_id", surf1.id_);
        }
      }
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

SurfaceXPlane::SurfaceXPlane(pugi::xml_node surf_node) : Surface(surf_node)
{
  read_coeffs(surf_node, id_, {&x0_});
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

SurfaceYPlane::SurfaceYPlane(pugi::xml_node surf_node) : Surface(surf_node)
{
  read_coeffs(surf_node, id_, {&y0_});
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

SurfaceZPlane::SurfaceZPlane(pugi::xml_node surf_node) : Surface(surf_node)
{
  read_coeffs(surf_node, id_, {&z0_});
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

SurfacePlane::SurfacePlane(pugi::xml_node surf_node) : Surface(surf_node)
{
  read_coeffs(surf_node, id_, {&A_, &B_, &C_, &D_});
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
  : Surface(surf_node)
{
  read_coeffs(surf_node, id_, {&y0_, &z0_, &radius_});
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
  : Surface(surf_node)
{
  read_coeffs(surf_node, id_, {&x0_, &z0_, &radius_});
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
  : Surface(surf_node)
{
  read_coeffs(surf_node, id_, {&x0_, &y0_, &radius_});
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

SurfaceSphere::SurfaceSphere(pugi::xml_node surf_node) : Surface(surf_node)
{
  read_coeffs(surf_node, id_, {&x0_, &y0_, &z0_, &radius_});
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

SurfaceXCone::SurfaceXCone(pugi::xml_node surf_node) : Surface(surf_node)
{
  read_coeffs(surf_node, id_, {&x0_, &y0_, &z0_, &radius_sq_});
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

SurfaceYCone::SurfaceYCone(pugi::xml_node surf_node) : Surface(surf_node)
{
  read_coeffs(surf_node, id_, {&x0_, &y0_, &z0_, &radius_sq_});
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

SurfaceZCone::SurfaceZCone(pugi::xml_node surf_node) : Surface(surf_node)
{
  read_coeffs(surf_node, id_, {&x0_, &y0_, &z0_, &radius_sq_});
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

SurfaceQuadric::SurfaceQuadric(pugi::xml_node surf_node) : Surface(surf_node)
{
  read_coeffs(
    surf_node, id_, {&A_, &B_, &C_, &D_, &E_, &F_, &G_, &H_, &J_, &K_});
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
// Torus helper functions
//==============================================================================

double torus_distance(double x1, double x2, double x3, double u1, double u2,
  double u3, double A, double B, double C, bool coincident)
{
  // Coefficients for equation: (c2 t^2 + c1 t + c0)^2 = c2' t^2 + c1' t + c0'
  double D = (C * C) / (B * B);
  double c2 = u1 * u1 + u2 * u2 + D * u3 * u3;
  double c1 = 2 * (u1 * x1 + u2 * x2 + D * u3 * x3);
  double c0 = x1 * x1 + x2 * x2 + D * x3 * x3 + A * A - C * C;
  double four_A2 = 4 * A * A;
  double c2p = four_A2 * (u1 * u1 + u2 * u2);
  double c1p = 2 * four_A2 * (u1 * x1 + u2 * x2);
  double c0p = four_A2 * (x1 * x1 + x2 * x2);

  // Coefficient for equation: a t^4 + b t^3 + c t^2 + d t + e = 0. If the point
  // is coincident, the 'e' coefficient should be zero. Explicitly setting it to
  // zero helps avoid numerical issues below with root finding.
  double coeff[5];
  coeff[0] = coincident ? 0.0 : c0 * c0 - c0p;
  coeff[1] = 2 * c0 * c1 - c1p;
  coeff[2] = c1 * c1 + 2 * c0 * c2 - c2p;
  coeff[3] = 2 * c1 * c2;
  coeff[4] = c2 * c2;

  std::complex<double> roots[4];
  oqs::quartic_solver(coeff, roots);

  // Find smallest positive, real root. In the case where the particle is
  // coincident with the surface, we are sure to have one root very close to
  // zero but possibly small and positive. A tolerance is set to discard that
  // zero.
  double distance = INFTY;
  double cutoff = coincident ? TORUS_TOL : 0.0;
  for (int i = 0; i < 4; ++i) {
    if (roots[i].imag() == 0) {
      double root = roots[i].real();
      if (root > cutoff && root < distance) {
        // Avoid roots corresponding to internal surfaces
        double s1 = x1 + u1 * root;
        double s2 = x2 + u2 * root;
        double s3 = x3 + u3 * root;
        double check = D * s3 * s3 + s1 * s1 + s2 * s2 + A * A - C * C;
        if (check >= 0) {
          distance = root;
        }
      }
    }
  }
  return distance;
}

//==============================================================================
// SurfaceXTorus implementation
//==============================================================================

SurfaceXTorus::SurfaceXTorus(pugi::xml_node surf_node) : Surface(surf_node)
{
  read_coeffs(surf_node, id_, {&x0_, &y0_, &z0_, &A_, &B_, &C_});
}

void SurfaceXTorus::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "x-torus", false);
  std::array<double, 6> coeffs {{x0_, y0_, z0_, A_, B_, C_}};
  write_dataset(group_id, "coefficients", coeffs);
}

double SurfaceXTorus::evaluate(Position r) const
{
  double x = r.x - x0_;
  double y = r.y - y0_;
  double z = r.z - z0_;
  return (x * x) / (B_ * B_) +
         std::pow(std::sqrt(y * y + z * z) - A_, 2) / (C_ * C_) - 1.;
}

double SurfaceXTorus::distance(Position r, Direction u, bool coincident) const
{
  double x = r.x - x0_;
  double y = r.y - y0_;
  double z = r.z - z0_;
  return torus_distance(y, z, x, u.y, u.z, u.x, A_, B_, C_, coincident);
}

Direction SurfaceXTorus::normal(Position r) const
{
  // reduce the expansion of the full form for torus
  double x = r.x - x0_;
  double y = r.y - y0_;
  double z = r.z - z0_;

  // f(x,y,z) = x^2/B^2 + (sqrt(y^2 + z^2) - A)^2/C^2 - 1
  // ∂f/∂x = 2x/B^2
  // ∂f/∂y = 2y(g - A)/(g*C^2) where g = sqrt(y^2 + z^2)
  // ∂f/∂z = 2z(g - A)/(g*C^2)
  // Multiplying by g*C^2*B^2 / 2 gives:
  double g = std::sqrt(y * y + z * z);
  double nx = C_ * C_ * g * x;
  double ny = y * (g - A_) * B_ * B_;
  double nz = z * (g - A_) * B_ * B_;
  Direction n(nx, ny, nz);
  return n / n.norm();
}

//==============================================================================
// SurfaceYTorus implementation
//==============================================================================

SurfaceYTorus::SurfaceYTorus(pugi::xml_node surf_node) : Surface(surf_node)
{
  read_coeffs(surf_node, id_, {&x0_, &y0_, &z0_, &A_, &B_, &C_});
}

void SurfaceYTorus::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "y-torus", false);
  std::array<double, 6> coeffs {{x0_, y0_, z0_, A_, B_, C_}};
  write_dataset(group_id, "coefficients", coeffs);
}

double SurfaceYTorus::evaluate(Position r) const
{
  double x = r.x - x0_;
  double y = r.y - y0_;
  double z = r.z - z0_;
  return (y * y) / (B_ * B_) +
         std::pow(std::sqrt(x * x + z * z) - A_, 2) / (C_ * C_) - 1.;
}

double SurfaceYTorus::distance(Position r, Direction u, bool coincident) const
{
  double x = r.x - x0_;
  double y = r.y - y0_;
  double z = r.z - z0_;
  return torus_distance(x, z, y, u.x, u.z, u.y, A_, B_, C_, coincident);
}

Direction SurfaceYTorus::normal(Position r) const
{
  // reduce the expansion of the full form for torus
  double x = r.x - x0_;
  double y = r.y - y0_;
  double z = r.z - z0_;

  // f(x,y,z) = y^2/B^2 + (sqrt(x^2 + z^2) - A)^2/C^2 - 1
  // ∂f/∂x = 2x(g - A)/(g*C^2) where g = sqrt(x^2 + z^2)
  // ∂f/∂y = 2y/B^2
  // ∂f/∂z = 2z(g - A)/(g*C^2)
  // Multiplying by g*C^2*B^2 / 2 gives:
  double g = std::sqrt(x * x + z * z);
  double nx = x * (g - A_) * B_ * B_;
  double ny = C_ * C_ * g * y;
  double nz = z * (g - A_) * B_ * B_;
  Direction n(nx, ny, nz);
  return n / n.norm();
}

//==============================================================================
// SurfaceZTorus implementation
//==============================================================================

SurfaceZTorus::SurfaceZTorus(pugi::xml_node surf_node) : Surface(surf_node)
{
  read_coeffs(surf_node, id_, {&x0_, &y0_, &z0_, &A_, &B_, &C_});
}

void SurfaceZTorus::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "z-torus", false);
  std::array<double, 6> coeffs {{x0_, y0_, z0_, A_, B_, C_}};
  write_dataset(group_id, "coefficients", coeffs);
}

double SurfaceZTorus::evaluate(Position r) const
{
  double x = r.x - x0_;
  double y = r.y - y0_;
  double z = r.z - z0_;
  return (z * z) / (B_ * B_) +
         std::pow(std::sqrt(x * x + y * y) - A_, 2) / (C_ * C_) - 1.;
}

double SurfaceZTorus::distance(Position r, Direction u, bool coincident) const
{
  double x = r.x - x0_;
  double y = r.y - y0_;
  double z = r.z - z0_;
  return torus_distance(x, y, z, u.x, u.y, u.z, A_, B_, C_, coincident);
}

Direction SurfaceZTorus::normal(Position r) const
{
  // reduce the expansion of the full form for torus
  double x = r.x - x0_;
  double y = r.y - y0_;
  double z = r.z - z0_;

  // f(x,y,z) = z^2/B^2 + (sqrt(x^2 + y^2) - A)^2/C^2 - 1
  // ∂f/∂x = 2x(g - A)/(g*C^2) where g = sqrt(x^2 + y^2)
  // ∂f/∂y = 2y(g - A)/(g*C^2)
  // ∂f/∂z = 2z/B^2
  // Multiplying by g*C^2*B^2 / 2 gives:
  double g = std::sqrt(x * x + y * y);
  double nx = x * (g - A_) * B_ * B_;
  double ny = y * (g - A_) * B_ * B_;
  double nz = C_ * C_ * g * z;
  Position n(nx, ny, nz);
  return n / n.norm();
}

//==============================================================================
// SurfaceRevolution implementation
//==============================================================================

SurfaceRevolution::SurfaceRevolution(pugi::xml_node surf_node)
  : Surface(surf_node)
{
  // Read all coefficients: x0, y0, z0, axis_code, n_points, r1, z1, r2, z2, ...
  auto coeffs = get_node_array<double>(surf_node, "coeffs");
  if (coeffs.size() < 7) {
    fatal_error(fmt::format(
      "Surface {} (revolution) requires at least 7 coefficients "
      "(origin, axis, n_points, and at least 2 profile points)",
      id_));
  }

  x0_ = coeffs[0];
  y0_ = coeffs[1];
  z0_ = coeffs[2];
  axis_ = static_cast<int>(coeffs[3]);
  int n_points = static_cast<int>(coeffs[4]);

  if (axis_ < 0 || axis_ > 2) {
    fatal_error(fmt::format(
      "Surface {} has invalid axis code {}. Must be 0 (x), 1 (y), or 2 (z).",
      id_, axis_));
  }

  if (n_points < 2) {
    fatal_error(fmt::format(
      "Surface {} requires at least 2 profile points, got {}.",
      id_, n_points));
  }

  size_t expected_size = 5 + 2 * n_points;
  if (coeffs.size() != expected_size) {
    fatal_error(fmt::format(
      "Surface {} expects {} coefficients for {} profile points, got {}.",
      id_, expected_size, n_points, coeffs.size()));
  }

  // Read profile points
  rz_.reserve(n_points);
  for (int i = 0; i < n_points; ++i) {
    double r = coeffs[5 + 2 * i];
    double z = coeffs[5 + 2 * i + 1];
    if (r < 0.0) {
      fatal_error(fmt::format(
        "Surface {} profile point {} has negative radius {}.",
        id_, i, r));
    }
    rz_.emplace_back(r, z);
  }
}

void SurfaceRevolution::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "revolution", false);

  // Write coefficients: x0, y0, z0, axis, n_points, r1, z1, r2, z2, ...
  vector<double> coeffs;
  coeffs.push_back(x0_);
  coeffs.push_back(y0_);
  coeffs.push_back(z0_);
  coeffs.push_back(static_cast<double>(axis_));
  coeffs.push_back(static_cast<double>(rz_.size()));
  for (const auto& point : rz_) {
    coeffs.push_back(point.first);
    coeffs.push_back(point.second);
  }
  write_dataset(group_id, "coefficients", coeffs);
}

void SurfaceRevolution::get_radial_and_axial(
  Position r, double& r_point, double& z_point) const
{
  double x = r.x - x0_;
  double y = r.y - y0_;
  double z = r.z - z0_;

  switch (axis_) {
  case 0: // x-axis
    r_point = std::sqrt(y * y + z * z);
    z_point = x;
    break;
  case 1: // y-axis
    r_point = std::sqrt(x * x + z * z);
    z_point = y;
    break;
  case 2: // z-axis
  default:
    r_point = std::sqrt(x * x + y * y);
    z_point = z;
    break;
  }
}

double SurfaceRevolution::find_segment_radius(double z_point) const
{
  int n = static_cast<int>(rz_.size());

  // Check if z is before the first point
  if (z_point <= rz_[0].second) {
    return rz_[0].first;
  }

  // Check if z is after the last point
  if (z_point >= rz_[n - 1].second) {
    return rz_[n - 1].first;
  }

  // Find the segment containing z_point
  for (int i = 0; i < n - 1; ++i) {
    double r1 = rz_[i].first;
    double z1 = rz_[i].second;
    double r2 = rz_[i + 1].first;
    double z2 = rz_[i + 1].second;

    if (z1 <= z_point && z_point <= z2) {
      // Linear interpolation
      if (std::abs(z2 - z1) < FP_COINCIDENT) {
        return std::min(r1, r2);
      }
      double t = (z_point - z1) / (z2 - z1);
      return r1 + t * (r2 - r1);
    }
  }

  return rz_[n - 1].first;
}

double SurfaceRevolution::evaluate(Position r) const
{
  double r_point, z_point;
  get_radial_and_axial(r, r_point, z_point);

  double r_surface = find_segment_radius(z_point);
  return r_point - r_surface;
}

// Helper function to compute intersection with a conical frustum segment
// This solves for the intersection of a ray with a cone frustum defined by
// (r1, z1) to (r2, z2) revolved around the specified axis.
double SurfaceRevolution::cone_distance(
  Position r, Direction u, bool coincident,
  double r1, double z1, double r2, double z2) const
{
  // Transform position to local coordinates
  double px = r.x - x0_;
  double py = r.y - y0_;
  double pz = r.z - z0_;

  // Get the radial and axial components based on axis
  double p_radial_sq, p_axial, u_radial_x, u_radial_y, u_axial;
  double p_radial_x, p_radial_y;

  switch (axis_) {
  case 0: // x-axis: radial is (y,z), axial is x
    p_radial_x = py;
    p_radial_y = pz;
    p_axial = px;
    u_radial_x = u.y;
    u_radial_y = u.z;
    u_axial = u.x;
    break;
  case 1: // y-axis: radial is (x,z), axial is y
    p_radial_x = px;
    p_radial_y = pz;
    p_axial = py;
    u_radial_x = u.x;
    u_radial_y = u.z;
    u_axial = u.y;
    break;
  case 2: // z-axis: radial is (x,y), axial is z
  default:
    p_radial_x = px;
    p_radial_y = py;
    p_axial = pz;
    u_radial_x = u.x;
    u_radial_y = u.y;
    u_axial = u.z;
    break;
  }

  p_radial_sq = p_radial_x * p_radial_x + p_radial_y * p_radial_y;

  // Check for degenerate segment (same z)
  double dz = z2 - z1;
  double dr = r2 - r1;

  if (std::abs(dz) < FP_COINCIDENT) {
    // This is essentially a disk or annulus - no surface intersection for cone
    return INFTY;
  }

  // For a cone frustum: sqrt(radial^2) = r1 + (axial - z1) * (r2 - r1) / (z2 - z1)
  // Squaring: radial^2 = (r1 + (axial - z1) * slope)^2 where slope = dr/dz
  // Expanding for ray (p + t*u):
  // (p_radial_x + t*u_radial_x)^2 + (p_radial_y + t*u_radial_y)^2 =
  //    (r1 + (p_axial + t*u_axial - z1) * slope)^2

  double slope = dr / dz;
  double base_r = r1 + (p_axial - z1) * slope;

  // Left side: |p_radial + t*u_radial|^2 = |p_radial|^2 + 2*t*(p_radial·u_radial) + t^2*|u_radial|^2
  double u_radial_sq = u_radial_x * u_radial_x + u_radial_y * u_radial_y;
  double p_dot_u_radial = p_radial_x * u_radial_x + p_radial_y * u_radial_y;

  // Right side: (base_r + t*u_axial*slope)^2 = base_r^2 + 2*base_r*t*u_axial*slope + t^2*(u_axial*slope)^2
  double u_axial_slope = u_axial * slope;

  // Quadratic: At^2 + Bt + C = 0
  double A = u_radial_sq - u_axial_slope * u_axial_slope;
  double B = 2.0 * (p_dot_u_radial - base_r * u_axial_slope);
  double C = p_radial_sq - base_r * base_r;

  if (coincident) {
    C = 0.0;
  }

  // Solve quadratic
  double quad = B * B - 4.0 * A * C;

  if (quad < 0.0) {
    return INFTY;
  }

  // Handle special case where A ≈ 0 (ray parallel to cone surface)
  if (std::abs(A) < FP_COINCIDENT) {
    if (std::abs(B) < FP_COINCIDENT) {
      return INFTY;
    }
    double t = -C / B;
    if (t > (coincident ? FP_COINCIDENT : 0.0)) {
      // Check if intersection is within segment bounds
      double z_hit = p_axial + t * u_axial;
      if (z_hit >= z1 && z_hit <= z2) {
        // Check we're on the positive side of the cone (not inside-out)
        double r_hit = base_r + t * u_axial_slope;
        if (r_hit >= 0.0) {
          return t;
        }
      }
    }
    return INFTY;
  }

  double sqrt_quad = std::sqrt(quad);
  double t1 = (-B - sqrt_quad) / (2.0 * A);
  double t2 = (-B + sqrt_quad) / (2.0 * A);

  double cutoff = coincident ? FP_COINCIDENT : 0.0;
  double distance = INFTY;

  // Check both roots
  for (double t : {t1, t2}) {
    if (t > cutoff && t < distance) {
      // Check if intersection is within segment bounds
      double z_hit = p_axial + t * u_axial;
      if (z_hit >= z1 - FP_COINCIDENT && z_hit <= z2 + FP_COINCIDENT) {
        // Check we're on the positive side of the cone
        double r_hit = base_r + t * u_axial_slope;
        if (r_hit >= -FP_COINCIDENT) {
          distance = t;
        }
      }
    }
  }

  return distance;
}

double SurfaceRevolution::distance(
  Position r, Direction u, bool coincident) const
{
  double min_distance = INFTY;
  int n = static_cast<int>(rz_.size());

  // Check intersection with each segment
  for (int i = 0; i < n - 1; ++i) {
    double r1 = rz_[i].first;
    double z1 = rz_[i].second;
    double r2 = rz_[i + 1].first;
    double z2 = rz_[i + 1].second;

    double d = cone_distance(r, u, coincident, r1, z1, r2, z2);
    if (d < min_distance) {
      min_distance = d;
    }
  }

  return min_distance;
}

Direction SurfaceRevolution::normal(Position r) const
{
  double r_point, z_point;
  get_radial_and_axial(r, r_point, z_point);

  // Find the segment containing this z position
  int n = static_cast<int>(rz_.size());
  int seg_idx = 0;

  if (z_point <= rz_[0].second) {
    seg_idx = 0;
  } else if (z_point >= rz_[n - 1].second) {
    seg_idx = n - 2;
  } else {
    for (int i = 0; i < n - 1; ++i) {
      if (z_point >= rz_[i].second && z_point <= rz_[i + 1].second) {
        seg_idx = i;
        break;
      }
    }
  }

  double r1 = rz_[seg_idx].first;
  double z1 = rz_[seg_idx].second;
  double r2 = rz_[seg_idx + 1].first;
  double z2 = rz_[seg_idx + 1].second;

  // For a cone: f = sqrt(radial^2) - (r1 + (axial - z1) * dr/dz) = 0
  // The gradient is:
  // ∂f/∂radial = radial / sqrt(radial^2) = radial / r_point (unit radial direction)
  // ∂f/∂axial = -dr/dz

  double dz = z2 - z1;
  double dr = r2 - r1;
  double slope = (std::abs(dz) > FP_COINCIDENT) ? (dr / dz) : 0.0;

  // Get position relative to origin
  double px = r.x - x0_;
  double py = r.y - y0_;
  double pz = r.z - z0_;

  // Compute the normal components
  double nx, ny, nz;

  if (r_point < FP_COINCIDENT) {
    // On the axis - normal points outward radially (arbitrary direction)
    // Use the axial gradient
    switch (axis_) {
    case 0:
      return Direction(-slope, 0.0, 0.0).normalize();
    case 1:
      return Direction(0.0, -slope, 0.0).normalize();
    case 2:
    default:
      return Direction(0.0, 0.0, -slope).normalize();
    }
  }

  // Radial unit vector components
  double radial_x, radial_y;

  switch (axis_) {
  case 0: // x-axis
    radial_x = py / r_point;
    radial_y = pz / r_point;
    // Normal = (radial_x, radial_y, 0) in (y, z, x) space
    // Mapping back: ny = radial_x, nz = radial_y, nx = -slope
    nx = -slope;
    ny = radial_x;
    nz = radial_y;
    break;
  case 1: // y-axis
    radial_x = px / r_point;
    radial_y = pz / r_point;
    // Normal in (x, z, y) space
    nx = radial_x;
    ny = -slope;
    nz = radial_y;
    break;
  case 2: // z-axis
  default:
    radial_x = px / r_point;
    radial_y = py / r_point;
    // Normal in (x, y, z) space
    nx = radial_x;
    ny = radial_y;
    nz = -slope;
    break;
  }

  Direction n(nx, ny, nz);
  double norm = n.norm();
  if (norm < FP_COINCIDENT) {
    // Degenerate case - return axial normal
    switch (axis_) {
    case 0:
      return Direction(1.0, 0.0, 0.0);
    case 1:
      return Direction(0.0, 1.0, 0.0);
    case 2:
    default:
      return Direction(0.0, 0.0, 1.0);
    }
  }

  return n / norm;
}

BoundingBox SurfaceRevolution::bounding_box(bool pos_side) const
{
  if (pos_side) {
    return {};
  }

  // Find the max radius and z extents
  double r_max = 0.0;
  double z_min = INFTY;
  double z_max = -INFTY;

  for (const auto& point : rz_) {
    r_max = std::max(r_max, point.first);
    z_min = std::min(z_min, point.second);
    z_max = std::max(z_max, point.second);
  }

  Position lower, upper;

  switch (axis_) {
  case 0: // x-axis
    lower = {x0_ + z_min, y0_ - r_max, z0_ - r_max};
    upper = {x0_ + z_max, y0_ + r_max, z0_ + r_max};
    break;
  case 1: // y-axis
    lower = {x0_ - r_max, y0_ + z_min, z0_ - r_max};
    upper = {x0_ + r_max, y0_ + z_max, z0_ + r_max};
    break;
  case 2: // z-axis
  default:
    lower = {x0_ - r_max, y0_ - r_max, z0_ + z_min};
    upper = {x0_ + r_max, y0_ + r_max, z0_ + z_max};
    break;
  }

  return {lower, upper};
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
  // periodic surfaces and their albedos.
  model::surfaces.reserve(n_surfaces);
  std::set<std::pair<int, int>> periodic_pairs;
  std::unordered_map<int, double> albedo_map;
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

      } else if (surf_type == "x-torus") {
        model::surfaces.push_back(std::make_unique<SurfaceXTorus>(surf_node));

      } else if (surf_type == "y-torus") {
        model::surfaces.push_back(std::make_unique<SurfaceYTorus>(surf_node));

      } else if (surf_type == "z-torus") {
        model::surfaces.push_back(std::make_unique<SurfaceZTorus>(surf_node));

      } else if (surf_type == "revolution") {
        model::surfaces.push_back(
          std::make_unique<SurfaceRevolution>(surf_node));

      } else {
        fatal_error(fmt::format("Invalid surface type, \"{}\"", surf_type));
      }

      // Check for a periodic surface
      if (check_for_node(surf_node, "boundary")) {
        std::string surf_bc = get_node_value(surf_node, "boundary", true, true);
        if (surf_bc == "periodic") {
          // Check for surface albedo. Skip sanity check as it is already done
          // in the Surface class's constructor.
          if (check_for_node(surf_node, "albedo")) {
            albedo_map[model::surfaces.back()->id_] =
              std::stod(get_node_value(surf_node, "albedo"));
          }
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

  // Assign the periodic boundary conditions with albedos
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
      surf1.bc_ = make_unique<TranslationalPeriodicBC>(i_surf, j_surf);
      surf2.bc_ = make_unique<TranslationalPeriodicBC>(i_surf, j_surf);
    } else {
      surf1.bc_ = make_unique<RotationalPeriodicBC>(i_surf, j_surf);
      surf2.bc_ = make_unique<RotationalPeriodicBC>(i_surf, j_surf);
    }

    // If albedo data is present in albedo map, set the boundary albedo.
    if (albedo_map.count(surf1.id_)) {
      surf1.bc_->set_albedo(albedo_map[surf1.id_]);
    }
    if (albedo_map.count(surf2.id_)) {
      surf2.bc_->set_albedo(albedo_map[surf2.id_]);
    }
  }
}

void free_memory_surfaces()
{
  model::surfaces.clear();
  model::surface_map.clear();
}

} // namespace openmc
