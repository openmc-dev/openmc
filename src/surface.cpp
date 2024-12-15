#include "openmc/surface.h"

#include <cmath>
#include <complex>
#include <initializer_list>
#include <set>
#include <utility>

#include <fmt/core.h>
#include <gsl/gsl-lite.hpp>

#include "openmc/array.h"
#include "openmc/container_util.h"
#include "openmc/error.h"
#include "openmc/external/quartic_solver.h"
#include "openmc/hdf5_interface.h"
#include "openmc/math_functions.h"
#include "openmc/random_lcg.h"
#include "openmc/search.h"
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

bool Surface::sense(Position r, Direction u, double t) const
{
  // Evaluate the surface equation at the particle's coordinates to determine
  // which side the particle is on.
  const double f = evaluate(r, t);

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
  Position r, Direction u, uint64_t* seed, GeometryState* p) const
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
      bc_->to_hdf5(surf_group);
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

  // Not moving?
  if (!(check_for_node(surf_node, "moving_velocities") ||
        check_for_node(surf_node, "moving_durations"))) {
    moving_ = false;
    return;
  }

  // Now, set the surface moving parameters
  moving_ = true;

  // Moving durations
  auto durations = get_node_array<double>(surf_node, "moving_durations");
  const int N_move = durations.size() + 1;

  // Moving time grids
  moving_time_grid_.resize(N_move + 1);
  moving_time_grid_[0] = 0.0;
  for (int n = 0; n < N_move - 1; n++) {
    moving_time_grid_[n + 1] = moving_time_grid_[n] + durations[n];
  }
  moving_time_grid_[N_move] = INFTY;

  // Moving velocities
  moving_velocities_.resize(N_move);
  std::string velocities_spec = get_node_value(surf_node, "moving_velocities");
  // Parse
  std::vector<double> numbers;
  for (int i = 0; i < velocities_spec.size();) {
    if (velocities_spec[i] == '-' || std::isdigit(velocities_spec[i])) {
      int j = i + 1;
      while (j < velocities_spec.size() && std::isdigit(velocities_spec[j])) {
        j++;
      }
      numbers.push_back(std::stod(velocities_spec.substr(i, j - i)));
      i = j;
    }
    i++;
  }
  // Assign to velocities
  for (int n = 0; n < N_move - 1; n++) {
    int idx = 3 * n;
    moving_velocities_[n][0] = numbers[idx];
    moving_velocities_[n][1] = numbers[idx + 1];
    moving_velocities_[n][2] = numbers[idx + 2];
  }
  moving_velocities_[N_move - 1] *= 0.0;

  // Moving translations
  moving_translations_.resize(N_move + 1);
  moving_translations_[0] *= 0.0;
  for (int n = 0; n < N_move - 1; n++) {
    moving_translations_[n + 1] =
      moving_translations_[n] + moving_velocities_[n] * durations[n];
  }
  moving_translations_[N_move] = moving_translations_[N_move - 1];
};

double CSGSurface::evaluate(Position r, double t) const
{
  if (!moving_) {
    return _evaluate(r);
  }
  // The surface moves

  // Get moving index
  int idx = lower_bound_index(
    moving_time_grid_.begin(), moving_time_grid_.end(), t);

  // Get moving translation, velocity, and starting time
  Position translation = moving_translations_[idx];
  double time_0 = moving_time_grid_[idx];
  Position velocity = moving_velocities_[idx];

  // Move the position relative to the surface movement
  double t_local = t - time_0;
  Position r_moved = r - (translation + velocity * t_local);
 
  // Evaluate the moved position
  return _evaluate(r_moved);
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

SurfaceXPlane::SurfaceXPlane(pugi::xml_node surf_node) : CSGSurface(surf_node)
{
  read_coeffs(surf_node, id_, {&x0_});
}

double SurfaceXPlane::_evaluate(Position r) const
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
  read_coeffs(surf_node, id_, {&y0_});
}

double SurfaceYPlane::_evaluate(Position r) const
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
  read_coeffs(surf_node, id_, {&z0_});
}

double SurfaceZPlane::_evaluate(Position r) const
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
  read_coeffs(surf_node, id_, {&A_, &B_, &C_, &D_});
}

double SurfacePlane::_evaluate(Position r) const
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
  read_coeffs(surf_node, id_, {&y0_, &z0_, &radius_});
}

double SurfaceXCylinder::_evaluate(Position r) const
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
  read_coeffs(surf_node, id_, {&x0_, &z0_, &radius_});
}

double SurfaceYCylinder::_evaluate(Position r) const
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
  read_coeffs(surf_node, id_, {&x0_, &y0_, &radius_});
}

double SurfaceZCylinder::_evaluate(Position r) const
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
  read_coeffs(surf_node, id_, {&x0_, &y0_, &z0_, &radius_});
}

double SurfaceSphere::_evaluate(Position r) const
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
  read_coeffs(surf_node, id_, {&x0_, &y0_, &z0_, &radius_sq_});
}

double SurfaceXCone::_evaluate(Position r) const
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
  read_coeffs(surf_node, id_, {&x0_, &y0_, &z0_, &radius_sq_});
}

double SurfaceYCone::_evaluate(Position r) const
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
  read_coeffs(surf_node, id_, {&x0_, &y0_, &z0_, &radius_sq_});
}

double SurfaceZCone::_evaluate(Position r) const
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
  read_coeffs(
    surf_node, id_, {&A_, &B_, &C_, &D_, &E_, &F_, &G_, &H_, &J_, &K_});
}

double SurfaceQuadric::_evaluate(Position r) const
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

SurfaceXTorus::SurfaceXTorus(pugi::xml_node surf_node) : CSGSurface(surf_node)
{
  read_coeffs(surf_node, id_, {&x0_, &y0_, &z0_, &A_, &B_, &C_});
}

void SurfaceXTorus::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "x-torus", false);
  std::array<double, 6> coeffs {{x0_, y0_, z0_, A_, B_, C_}};
  write_dataset(group_id, "coefficients", coeffs);
}

double SurfaceXTorus::_evaluate(Position r) const
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

SurfaceYTorus::SurfaceYTorus(pugi::xml_node surf_node) : CSGSurface(surf_node)
{
  read_coeffs(surf_node, id_, {&x0_, &y0_, &z0_, &A_, &B_, &C_});
}

void SurfaceYTorus::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "y-torus", false);
  std::array<double, 6> coeffs {{x0_, y0_, z0_, A_, B_, C_}};
  write_dataset(group_id, "coefficients", coeffs);
}

double SurfaceYTorus::_evaluate(Position r) const
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

SurfaceZTorus::SurfaceZTorus(pugi::xml_node surf_node) : CSGSurface(surf_node)
{
  read_coeffs(surf_node, id_, {&x0_, &y0_, &z0_, &A_, &B_, &C_});
}

void SurfaceZTorus::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "z-torus", false);
  std::array<double, 6> coeffs {{x0_, y0_, z0_, A_, B_, C_}};
  write_dataset(group_id, "coefficients", coeffs);
}

double SurfaceZTorus::_evaluate(Position r) const
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
