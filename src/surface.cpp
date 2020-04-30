#include "openmc/surface.h"

#include <array>
#include <cmath>
#include <utility>

#include <fmt/core.h>

#include "openmc/error.h"
#include "openmc/dagmc.h"
#include "openmc/hdf5_interface.h"
#include "openmc/settings.h"
#include "openmc/string_utils.h"
#include "openmc/xml_interface.h"
#include "openmc/random_lcg.h"
#include "openmc/math_functions.h"

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace model {
  //std::vector<std::unique_ptr<Surface>> surfaces;
  std::vector<Surface> surfaces;
  Surface* device_surfaces;
  std::unordered_map<int, int> surface_map;
} // namespace model

//==============================================================================
// Helper functions for reading the "coeffs" node of an XML surface element
//==============================================================================

void read_coeffs(pugi::xml_node surf_node, int surf_id, double &c1)
{
  // Check the given number of coefficients.
  std::string coeffs = get_node_value(surf_node, "coeffs");
  int n_words = word_count(coeffs);
  if (n_words != 1) {
    fatal_error(fmt::format("Surface {} expects 1 coeff but was given {}",
      surf_id, n_words));
  }

  // Parse the coefficients.
  int stat = sscanf(coeffs.c_str(), "%lf", &c1);
  if (stat != 1) {
    fatal_error("Something went wrong reading surface coeffs");
  }
}

void read_coeffs(pugi::xml_node surf_node, int surf_id, double &c1, double &c2,
                 double &c3)
{
  // Check the given number of coefficients.
  std::string coeffs = get_node_value(surf_node, "coeffs");
  int n_words = word_count(coeffs);
  if (n_words != 3) {
    fatal_error(fmt::format("Surface {} expects 3 coeffs but was given {}",
      surf_id, n_words));
  }

  // Parse the coefficients.
  int stat = sscanf(coeffs.c_str(), "%lf %lf %lf", &c1, &c2, &c3);
  if (stat != 3) {
    fatal_error("Something went wrong reading surface coeffs");
  }
}

void read_coeffs(pugi::xml_node surf_node, int surf_id, double &c1, double &c2,
                 double &c3, double &c4)
{
  // Check the given number of coefficients.
  std::string coeffs = get_node_value(surf_node, "coeffs");
  int n_words = word_count(coeffs);
  if (n_words != 4) {
    fatal_error(fmt::format("Surface {} expects 4 coeffs but was given ",
      surf_id, n_words));
  }

  // Parse the coefficients.
  int stat = sscanf(coeffs.c_str(), "%lf %lf %lf %lf", &c1, &c2, &c3, &c4);
  if (stat != 4) {
    fatal_error("Something went wrong reading surface coeffs");
  }
}

void read_coeffs(pugi::xml_node surf_node, int surf_id, double &c1, double &c2,
                 double &c3, double &c4, double &c5, double &c6, double &c7,
                 double &c8, double &c9, double &c10)
{
  // Check the given number of coefficients.
  std::string coeffs = get_node_value(surf_node, "coeffs");
  int n_words = word_count(coeffs);
  if (n_words != 10) {
    fatal_error(fmt::format("Surface {} expects 10 coeffs but was given {}",
      surf_id, n_words));
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

Surface::Surface(pugi::xml_node surf_node, Surface::SurfaceType type) : type_(type)
{
  if (check_for_node(surf_node, "id")) {
    id_ = std::stoi(get_node_value(surf_node, "id"));
  } else {
    fatal_error("Must specify id of surface in geometry XML file.");
  }

  if (check_for_node(surf_node, "name")) {
    name_ = get_node_value(surf_node, "name", false);
  }

  if (check_for_node(surf_node, "boundary")) {
    std::string surf_bc = get_node_value(surf_node, "boundary", true, true);

    if (surf_bc == "transmission" || surf_bc == "transmit" ||surf_bc.empty()) {
      bc_ = BoundaryType::TRANSMIT;

    } else if (surf_bc == "vacuum") {
      bc_ = BoundaryType::VACUUM;

    } else if (surf_bc == "reflective" || surf_bc == "reflect"
               || surf_bc == "reflecting") {
      bc_ = BoundaryType::REFLECT;
    } else if (surf_bc == "white") {
      bc_ = BoundaryType::WHITE;
    } else if (surf_bc == "periodic") {
      bc_ = BoundaryType::PERIODIC;
    } else {
      fatal_error(fmt::format("Unknown boundary condition \"{}\" specified "
        "on surface {}", surf_bc, id_));
    }

  } else {
    bc_ = BoundaryType::TRANSMIT;
  }
  if (check_for_node(surf_node, "periodic_surface_id")) {
    i_periodic_ = std::stoi(get_node_value(surf_node, "periodic_surface_id"));
  }

  // Type specific initializations
  switch(type_){
    case SurfaceType::SurfaceXPlane:    read_coeffs(surf_node, id_, x0_); break;
    case SurfaceType::SurfaceYPlane:    read_coeffs(surf_node, id_, y0_); break;
    case SurfaceType::SurfaceZPlane:    read_coeffs(surf_node, id_, z0_); break;
    case SurfaceType::SurfacePlane:     read_coeffs(surf_node, id_, A_, B_, C_, D_); break;
    case SurfaceType::SurfaceXCylinder: read_coeffs(surf_node, id_, y0_, z0_, radius_); break;
    case SurfaceType::SurfaceYCylinder: read_coeffs(surf_node, id_, x0_, z0_, radius_); break;
    case SurfaceType::SurfaceZCylinder: read_coeffs(surf_node, id_, x0_, y0_, radius_); break;
    case SurfaceType::SurfaceSphere:    read_coeffs(surf_node, id_, x0_, y0_, z0_, radius_); break;
    case SurfaceType::SurfaceXCone:     read_coeffs(surf_node, id_, x0_, y0_, z0_, radius_sq_); break;
    case SurfaceType::SurfaceYCone:     read_coeffs(surf_node, id_, x0_, y0_, z0_, radius_sq_); break;
    case SurfaceType::SurfaceZCone:     read_coeffs(surf_node, id_, x0_, y0_, z0_, radius_sq_); break;
    case SurfaceType::SurfaceQuadric:   read_coeffs(surf_node, id_, A_, B_, C_, D_, E_, F_, G_, H_, J_, K_); break;
  }
}

bool
Surface::sense(Position r, Direction u) const
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

Direction
Surface::reflect(Position r, Direction u, Particle* p) const
{
  // Determine projection of direction onto normal and squared magnitude of
  // normal.
  Direction n = normal(r);
  const double projection = n.dot(u);
  const double magnitude = n.dot(n);

  // Reflect direction according to normal.
  return u -= (2.0 * projection / magnitude) * n;
}

Direction
Surface::diffuse_reflect(Position r, Direction u, uint64_t* seed) const
{
  // Diffuse reflect direction according to the normal.
  // cosine distribution

  Direction n = this->normal(r);
  n /= n.norm();
  const double projection = n.dot(u);

  // sample from inverse function, u=sqrt(rand) since p(u)=2u, so F(u)=u^2
  const double mu = (projection>=0.0) ?
                  -std::sqrt(prn(seed)) : std::sqrt(prn(seed));

  // sample azimuthal distribution uniformly
  u = rotate_angle(n, mu, nullptr, seed);

  // normalize the direction
  return u/u.norm();
}


void
Surface::to_hdf5(hid_t group_id) const
{
  std::string group_name {"surface "};
  group_name += std::to_string(id_);

  hid_t surf_group = create_group(group_id, group_name);

  switch(bc_) {
    case BoundaryType::TRANSMIT :
      write_string(surf_group, "boundary_type", "transmission", false);
      break;
    case BoundaryType::VACUUM :
      write_string(surf_group, "boundary_type", "vacuum", false);
      break;
    case BoundaryType::REFLECT :
      write_string(surf_group, "boundary_type", "reflective", false);
      break;
    case BoundaryType::WHITE :
      write_string(surf_group, "boundary_type", "white", false);
      break;
    case BoundaryType::PERIODIC :
      write_string(surf_group, "boundary_type", "periodic", false);
      break;
  }

  if (!name_.empty()) {
    write_string(surf_group, "name", name_, false);
  }

  to_hdf5_inner(surf_group);

  close_group(surf_group);
}


//==============================================================================
// DAGSurface implementation
//==============================================================================
/*
#ifdef DAGMC
DAGSurface::DAGSurface() : Surface{} {} // empty constructor

double DAGSurface::evaluate(Position r) const
{
  return 0.0;
}

double
DAGSurface::distance(Position r, Direction u, bool coincident) const
{
  moab::ErrorCode rval;
  moab::EntityHandle surf = dagmc_ptr_->entity_by_index(2, dag_index_);
  moab::EntityHandle hit_surf;
  double dist;
  double pnt[3] = {r.x, r.y, r.z};
  double dir[3] = {u.x, u.y, u.z};
  rval = dagmc_ptr_->ray_fire(surf, pnt, dir, hit_surf, dist, NULL, 0, 0);
  MB_CHK_ERR_CONT(rval);
  if (dist < 0.0) dist = INFTY;
  return dist;
}

Direction DAGSurface::normal(Position r) const
{
  moab::ErrorCode rval;
  moab::EntityHandle surf = dagmc_ptr_->entity_by_index(2, dag_index_);
  double pnt[3] = {r.x, r.y, r.z};
  double dir[3];
  rval = dagmc_ptr_->get_angle(surf, pnt, dir);
  MB_CHK_ERR_CONT(rval);
  return dir;
}

Direction DAGSurface::reflect(Position r, Direction u, Particle* p) const
{
  p->history_.reset_to_last_intersection();
  p->last_dir_ = Surface::reflect(r, u, p);
  return p->last_dir_;
}

void DAGSurface::to_hdf5(hid_t group_id) const {}

#endif
*/
//==============================================================================
// PeriodicSurface implementation
//==============================================================================

/*
PeriodicSurface::PeriodicSurface(pugi::xml_node surf_node)
  : Surface {surf_node}
{
  if (check_for_node(surf_node, "periodic_surface_id")) {
    i_periodic_ = std::stoi(get_node_value(surf_node, "periodic_surface_id"));
  }
}
*/

//==============================================================================
// Generic functions for x-, y-, and z-, planes.
//==============================================================================

// The template parameter indicates the axis normal to the plane.
template<int i> double
axis_aligned_plane_distance(Position r, Direction u, bool coincident, double offset)
{
  const double f = offset - r[i];
  if (coincident || std::abs(f) < FP_COINCIDENT || u[i] == 0.0) return INFTY;
  const double d = f / u[i];
  if (d < 0.0) return INFTY;
  return d;
}

//==============================================================================
// SurfaceXPlane implementation
//==============================================================================

double Surface::SurfaceXPlane_evaluate(Position r) const
{
  return r.x - x0_;
}

double Surface::SurfaceXPlane_distance(Position r, Direction u, bool coincident) const
{
  return axis_aligned_plane_distance<0>(r, u, coincident, x0_);
}

Direction Surface::SurfaceXPlane_normal(Position r) const
{
  return {1., 0., 0.};
}

void Surface::SurfaceXPlane_to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "x-plane", false);
  std::array<double, 1> coeffs {{x0_}};
  write_dataset(group_id, "coefficients", coeffs);
}

bool Surface::SurfaceXPlane_periodic_translate(const Surface* other,
                                       Position& r,  Direction& u) const
{
  Direction other_n = other->normal(r);
  if (other_n.x == 1 && other_n.y == 0 && other_n.z == 0) {
    r.x = x0_;
    return false;
  } else {
    // Assume the partner is an YPlane (the only supported partner).  Use the
    // evaluate function to find y0, then adjust position/Direction for
    // rotational symmetry.
    double y0_ = -other->evaluate({0., 0., 0.});
    r.y = r.x - x0_ + y0_;
    r.x = x0_;

    double ux = u.x;
    u.x = -u.y;
    u.y = ux;

    return true;
  }
}

BoundingBox
Surface::SurfaceXPlane_bounding_box(bool pos_side) const
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

double Surface::SurfaceYPlane_evaluate(Position r) const
{
  return r.y - y0_;
}

double Surface::SurfaceYPlane_distance(Position r, Direction u, bool coincident) const
{
  return axis_aligned_plane_distance<1>(r, u, coincident, y0_);
}

Direction Surface::SurfaceYPlane_normal(Position r) const
{
  return {0., 1., 0.};
}

void Surface::SurfaceYPlane_to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "y-plane", false);
  std::array<double, 1> coeffs {{y0_}};
  write_dataset(group_id, "coefficients", coeffs);
}

bool Surface::SurfaceYPlane_periodic_translate(const Surface* other,
                                       Position& r, Direction& u) const
{
  Direction other_n = other->normal(r);
  if (other_n.x == 0 && other_n.y == 1 && other_n.z == 0) {
    // The periodic partner is also aligned along y.  Just change the y coord.
    r.y = y0_;
    return false;
  } else {
    // Assume the partner is an XPlane (the only supported partner).  Use the
    // evaluate function to find x0, then adjust position/Direction for rotational
    // symmetry.
    double x0_ = -other->evaluate({0., 0., 0.});
    r.x = r.y - y0_ + x0_;
    r.y = y0_;

    double ux = u.x;
    u.x = u.y;
    u.y = -ux;

    return true;
  }
}

BoundingBox
Surface::SurfaceYPlane_bounding_box(bool pos_side) const
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

double Surface::SurfaceZPlane_evaluate(Position r) const
{
  return r.z - z0_;
}

double Surface::SurfaceZPlane_distance(Position r, Direction u, bool coincident) const
{
  return axis_aligned_plane_distance<2>(r, u, coincident, z0_);
}

Direction Surface::SurfaceZPlane_normal(Position r) const
{
  return {0., 0., 1.};
}

void Surface::SurfaceZPlane_to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "z-plane", false);
  std::array<double, 1> coeffs {{z0_}};
  write_dataset(group_id, "coefficients", coeffs);
}

bool Surface::SurfaceZPlane_periodic_translate(const Surface* other,
                                       Position& r, Direction& u) const
{
  // Assume the other plane is aligned along z.  Just change the z coord.
  r.z = z0_;
  return false;
}

BoundingBox
Surface::SurfaceZPlane_bounding_box(bool pos_side) const
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

double
Surface::SurfacePlane_evaluate(Position r) const
{
  return A_*r.x + B_*r.y + C_*r.z - D_;
}

double
Surface::SurfacePlane_distance(Position r, Direction u, bool coincident) const
{
  const double f = A_*r.x + B_*r.y + C_*r.z - D_;
  const double projection = A_*u.x + B_*u.y + C_*u.z;
  if (coincident || std::abs(f) < FP_COINCIDENT || projection == 0.0) {
    return INFTY;
  } else {
    const double d = -f / projection;
    if (d < 0.0) return INFTY;
    return d;
  }
}

Direction
Surface::SurfacePlane_normal(Position r) const
{
  return {A_, B_, C_};
}

void Surface::SurfacePlane_to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "plane", false);
  std::array<double, 4> coeffs {{A_, B_, C_, D_}};
  write_dataset(group_id, "coefficients", coeffs);
}

bool Surface::SurfacePlane_periodic_translate(const Surface* other, Position& r,
                                      Direction& u) const
{
  // This function assumes the other plane shares this plane's normal direction.

  // Determine the distance to intersection.
  double d = evaluate(r) / (A_*A_ + B_*B_ + C_*C_);

  // Move the particle that distance along the normal vector.
  r.x -= d * A_;
  r.y -= d * B_;
  r.z -= d * C_;

  return false;
}

//==============================================================================
// Generic functions for x-, y-, and z-, cylinders
//==============================================================================

// The template parameters indicate the axes perpendicular to the axis of the
// cylinder.  offset1 and offset2 should correspond with i1 and i2,
// respectively.
template<int i1, int i2> double
axis_aligned_cylinder_evaluate(Position r, double offset1,
                               double offset2, double radius)
{
  const double r1 = r[i1] - offset1;
  const double r2 = r[i2] - offset2;
  return r1*r1 + r2*r2 - radius*radius;
}

// The first template parameter indicates which axis the cylinder is aligned to.
// The other two parameters indicate the other two axes.  offset1 and offset2
// should correspond with i2 and i3, respectively.
template<int i1, int i2, int i3> double
axis_aligned_cylinder_distance(Position r, Direction u,
     bool coincident, double offset1, double offset2, double radius)
{
  const double a = 1.0 - u[i1]*u[i1];  // u^2 + v^2
  if (a == 0.0) return INFTY;

  const double r2 = r[i2] - offset1;
  const double r3 = r[i3] - offset2;
  const double k = r2 * u[i2] + r3 * u[i3];
  const double c = r2*r2 + r3*r3 - radius*radius;
  const double quad = k*k - a*c;

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
    if (d < 0.0) return INFTY;
    return d;
  }
}

// The first template parameter indicates which axis the cylinder is aligned to.
// The other two parameters indicate the other two axes.  offset1 and offset2
// should correspond with i2 and i3, respectively.
template<int i1, int i2, int i3> Direction
axis_aligned_cylinder_normal(Position r, double offset1, double offset2)
{
  Direction u;
  u[i2] = 2.0 * (r[i2] - offset1);
  u[i3] = 2.0 * (r[i3] - offset2);
  u[i1] = 0.0;
  return u;
}

//==============================================================================
// SurfaceXCylinder implementation
//==============================================================================

double Surface::SurfaceXCylinder_evaluate(Position r) const
{
  return axis_aligned_cylinder_evaluate<1, 2>(r, y0_, z0_, radius_);
}

double Surface::SurfaceXCylinder_distance(Position r, Direction u, bool coincident) const
{
  return axis_aligned_cylinder_distance<0, 1, 2>(r, u, coincident, y0_, z0_,
                                                 radius_);
}

Direction Surface::SurfaceXCylinder_normal(Position r) const
{
  return axis_aligned_cylinder_normal<0, 1, 2>(r, y0_, z0_);
}

void Surface::SurfaceXCylinder_to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "x-cylinder", false);
  std::array<double, 3> coeffs {{y0_, z0_, radius_}};
  write_dataset(group_id, "coefficients", coeffs);
}

BoundingBox Surface::SurfaceXCylinder_bounding_box(bool pos_side) const {
  if (!pos_side) {
    return {-INFTY, INFTY, y0_ - radius_, y0_ + radius_, z0_ - radius_, z0_ + radius_};
  } else {
    return {};
  }
}
//==============================================================================
// SurfaceYCylinder implementation
//==============================================================================

double Surface::SurfaceYCylinder_evaluate(Position r) const
{
  return axis_aligned_cylinder_evaluate<0, 2>(r, x0_, z0_, radius_);
}

double Surface::SurfaceYCylinder_distance(Position r, Direction u, bool coincident) const
{
  return axis_aligned_cylinder_distance<1, 0, 2>(r, u, coincident, x0_, z0_,
                                                 radius_);
}

Direction Surface::SurfaceYCylinder_normal(Position r) const
{
  return axis_aligned_cylinder_normal<1, 0, 2>(r, x0_, z0_);
}

void Surface::SurfaceYCylinder_to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "y-cylinder", false);
  std::array<double, 3> coeffs {{x0_, z0_, radius_}};
  write_dataset(group_id, "coefficients", coeffs);
}

BoundingBox Surface::SurfaceYCylinder_bounding_box(bool pos_side) const {
  if (!pos_side) {
    return {x0_ - radius_, x0_ + radius_, -INFTY, INFTY, z0_ - radius_, z0_ + radius_};
  } else {
    return {};
  }
}

//==============================================================================
// SurfaceZCylinder implementation
//==============================================================================

double Surface::SurfaceZCylinder_evaluate(Position r) const
{
  return axis_aligned_cylinder_evaluate<0, 1>(r, x0_, y0_, radius_);
}

double Surface::SurfaceZCylinder_distance(Position r, Direction u, bool coincident) const
{
  return axis_aligned_cylinder_distance<2, 0, 1>(r, u, coincident, x0_, y0_,
                                                 radius_);
}

Direction Surface::SurfaceZCylinder_normal(Position r) const
{
  return axis_aligned_cylinder_normal<2, 0, 1>(r, x0_, y0_);
}

void Surface::SurfaceZCylinder_to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "z-cylinder", false);
  std::array<double, 3> coeffs {{x0_, y0_, radius_}};
  write_dataset(group_id, "coefficients", coeffs);
}

BoundingBox Surface::SurfaceZCylinder_bounding_box(bool pos_side) const {
  if (!pos_side) {
    return {x0_ - radius_, x0_ + radius_, y0_ - radius_, y0_ + radius_, -INFTY, INFTY};
  } else {
    return {};
  }
}


//==============================================================================
// SurfaceSphere implementation
//==============================================================================

double Surface::SurfaceSphere_evaluate(Position r) const
{
  const double x = r.x - x0_;
  const double y = r.y - y0_;
  const double z = r.z - z0_;
  return x*x + y*y + z*z - radius_*radius_;
}

double Surface::SurfaceSphere_distance(Position r, Direction u, bool coincident) const
{
  const double x = r.x - x0_;
  const double y = r.y - y0_;
  const double z = r.z - z0_;
  const double k = x*u.x + y*u.y + z*u.z;
  const double c = x*x + y*y + z*z - radius_*radius_;
  const double quad = k*k - c;

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
    if (d < 0.0) return INFTY;
    return d;
  }
}

Direction Surface::SurfaceSphere_normal(Position r) const
{
  return {2.0*(r.x - x0_), 2.0*(r.y - y0_), 2.0*(r.z - z0_)};
}

void Surface::SurfaceSphere_to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "sphere", false);
  std::array<double, 4> coeffs {{x0_, y0_, z0_, radius_}};
  write_dataset(group_id, "coefficients", coeffs);
}

BoundingBox Surface::SurfaceSphere_bounding_box(bool pos_side) const {
  if (!pos_side) {
    return {x0_ - radius_, x0_ + radius_,
            y0_ - radius_, y0_ + radius_,
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
template<int i1, int i2, int i3> double
axis_aligned_cone_evaluate(Position r, double offset1,
                           double offset2, double offset3, double radius_sq)
{
  const double r1 = r[i1] - offset1;
  const double r2 = r[i2] - offset2;
  const double r3 = r[i3] - offset3;
  return r2*r2 + r3*r3 - radius_sq*r1*r1;
}

// The first template parameter indicates which axis the cone is aligned to.
// The other two parameters indicate the other two axes.  offset1, offset2,
// and offset3 should correspond with i1, i2, and i3, respectively.
template<int i1, int i2, int i3> double
axis_aligned_cone_distance(Position r, Direction u,
     bool coincident, double offset1, double offset2, double offset3,
     double radius_sq)
{
  const double r1 = r[i1] - offset1;
  const double r2 = r[i2] - offset2;
  const double r3 = r[i3] - offset3;
  const double a = u[i2]*u[i2] + u[i3]*u[i3]
                   - radius_sq*u[i1]*u[i1];
  const double k = r2*u[i2] + r3*u[i3] - radius_sq*r1*u[i1];
  const double c = r2*r2 + r3*r3 - radius_sq*r1*r1;
  double quad = k*k - a*c;

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
      if (b > 0.0) d = b;
    } else {
      if (b > 0.0) {
        if (b < d) d = b;
      }
    }
  }

  // If the distance was negative, set boundary distance to infinity.
  if (d <= 0.0) return INFTY;
  return d;
}

// The first template parameter indicates which axis the cone is aligned to.
// The other two parameters indicate the other two axes.  offset1, offset2,
// and offset3 should correspond with i1, i2, and i3, respectively.
template<int i1, int i2, int i3> Direction
axis_aligned_cone_normal(Position r, double offset1, double offset2,
                         double offset3, double radius_sq)
{
  Direction u;
  u[i1] = -2.0 * radius_sq * (r[i1] - offset1);
  u[i2] = 2.0 * (r[i2] - offset2);
  u[i3] = 2.0 * (r[i3] - offset3);
  return u;
}

//==============================================================================
// SurfaceXCone implementation
//==============================================================================

double Surface::SurfaceXCone_evaluate(Position r) const
{
  return axis_aligned_cone_evaluate<0, 1, 2>(r, x0_, y0_, z0_, radius_sq_);
}

double Surface::SurfaceXCone_distance(Position r, Direction u, bool coincident) const
{
  return axis_aligned_cone_distance<0, 1, 2>(r, u, coincident, x0_, y0_, z0_,
                                             radius_sq_);
}

Direction Surface::SurfaceXCone_normal(Position r) const
{
  return axis_aligned_cone_normal<0, 1, 2>(r, x0_, y0_, z0_, radius_sq_);
}

void Surface::SurfaceXCone_to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "x-cone", false);
  std::array<double, 4> coeffs {{x0_, y0_, z0_, radius_sq_}};
  write_dataset(group_id, "coefficients", coeffs);
}

//==============================================================================
// SurfaceYCone implementation
//==============================================================================

double Surface::SurfaceYCone_evaluate(Position r) const
{
  return axis_aligned_cone_evaluate<1, 0, 2>(r, y0_, x0_, z0_, radius_sq_);
}

double Surface::SurfaceYCone_distance(Position r, Direction u, bool coincident) const
{
  return axis_aligned_cone_distance<1, 0, 2>(r, u, coincident, y0_, x0_, z0_,
                                             radius_sq_);
}

Direction Surface::SurfaceYCone_normal(Position r) const
{
  return axis_aligned_cone_normal<1, 0, 2>(r, y0_, x0_, z0_, radius_sq_);
}

void Surface::SurfaceYCone_to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "y-cone", false);
  std::array<double, 4> coeffs {{x0_, y0_, z0_, radius_sq_}};
  write_dataset(group_id, "coefficients", coeffs);
}

//==============================================================================
// SurfaceZCone implementation
//==============================================================================

double Surface::SurfaceZCone_evaluate(Position r) const
{
  return axis_aligned_cone_evaluate<2, 0, 1>(r, z0_, x0_, y0_, radius_sq_);
}

double Surface::SurfaceZCone_distance(Position r, Direction u, bool coincident) const
{
  return axis_aligned_cone_distance<2, 0, 1>(r, u, coincident, z0_, x0_, y0_,
                                             radius_sq_);
}

Direction Surface::SurfaceZCone_normal(Position r) const
{
  return axis_aligned_cone_normal<2, 0, 1>(r, z0_, x0_, y0_, radius_sq_);
}

void Surface::SurfaceZCone_to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "z-cone", false);
  std::array<double, 4> coeffs {{x0_, y0_, z0_, radius_sq_}};
  write_dataset(group_id, "coefficients", coeffs);
}

//==============================================================================
// SurfaceQuadric implementation
//==============================================================================

double
Surface::SurfaceQuadric_evaluate(Position r) const
{
  const double x = r.x;
  const double y = r.y;
  const double z = r.z;
  return x*(A_*x + D_*y + G_) +
         y*(B_*y + E_*z + H_) +
         z*(C_*z + F_*x + J_) + K_;
}

double
Surface::SurfaceQuadric_distance(Position r, Direction ang, bool coincident) const
{
  const double &x = r.x;
  const double &y = r.y;
  const double &z = r.z;
  const double &u = ang.x;
  const double &v = ang.y;
  const double &w = ang.z;

  const double a = A_*u*u + B_*v*v + C_*w*w + D_*u*v + E_*v*w + F_*u*w;
  const double k = A_*u*x + B_*v*y + C_*w*z + 0.5*(D_*(u*y + v*x)
                   + E_*(v*z + w*y) + F_*(w*x + u*z) + G_*u + H_*v + J_*w);
  const double c = A_*x*x + B_*y*y + C_*z*z + D_*x*y + E_*y*z +  F_*x*z + G_*x
                   + H_*y + J_*z + K_;
  double quad = k*k - a*c;

  double d;

  if (quad < 0.0) {
    // No intersection with surface.
    return INFTY;

  } else if (coincident || std::abs(c) < FP_COINCIDENT) {
    // Particle is on the surface, thus one distance is positive/negative and
    // the other is zero. The sign of k determines which distance is zero and
    // which is not.
    if (k >= 0.0) {
      d = (-k - sqrt(quad)) / a;
    } else {
      d = (-k + sqrt(quad)) / a;
    }

  } else {
    // Calculate both solutions to the quadratic.
    quad = sqrt(quad);
    d = (-k - quad) / a;
    double b = (-k + quad) / a;

    // Determine the smallest positive solution.
    if (d < 0.0) {
      if (b > 0.0) d = b;
    } else {
      if (b > 0.0) {
        if (b < d) d = b;
      }
    }
  }

  // If the distance was negative, set boundary distance to infinity.
  if (d <= 0.0) return INFTY;
  return d;
}

Direction
Surface::SurfaceQuadric_normal(Position r) const
{
  const double &x = r.x;
  const double &y = r.y;
  const double &z = r.z;
  return {2.0*A_*x + D_*y + F_*z + G_,
          2.0*B_*y + D_*x + E_*z + H_,
          2.0*C_*z + E_*y + F_*x + J_};
}

void Surface::SurfaceQuadric_to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "quadric", false);
  std::array<double, 10> coeffs {{A_, B_, C_, D_, E_, F_, G_, H_, J_, K_}};
  write_dataset(group_id, "coefficients", coeffs);
}

//==============================================================================
// Dispatchers
//==============================================================================
double Surface::evaluate(Position r) const
{
  switch(type_){
    case SurfaceType::SurfaceXPlane:    return SurfaceXPlane_evaluate(r);    break;
    case SurfaceType::SurfaceYPlane:    return SurfaceYPlane_evaluate(r);    break;
    case SurfaceType::SurfaceZPlane:    return SurfaceZPlane_evaluate(r);    break;
    case SurfaceType::SurfacePlane:     return SurfacePlane_evaluate(r);     break;
    case SurfaceType::SurfaceXCylinder: return SurfaceXCylinder_evaluate(r); break;
    case SurfaceType::SurfaceYCylinder: return SurfaceYCylinder_evaluate(r); break;
    case SurfaceType::SurfaceZCylinder: return SurfaceZCylinder_evaluate(r); break;
    case SurfaceType::SurfaceSphere:    return SurfaceSphere_evaluate(r);    break;
    case SurfaceType::SurfaceXCone:     return SurfaceXCone_evaluate(r);     break;
    case SurfaceType::SurfaceYCone:     return SurfaceYCone_evaluate(r);     break;
    case SurfaceType::SurfaceZCone:     return SurfaceZCone_evaluate(r);     break;
    case SurfaceType::SurfaceQuadric:   return SurfaceQuadric_evaluate(r);   break;
  }
}

double Surface::distance(Position r, Direction u, bool coincident) const
{
  switch(type_){
    case SurfaceType::SurfaceXPlane:    return SurfaceXPlane_distance(r, u, coincident);    break;
    case SurfaceType::SurfaceYPlane:    return SurfaceYPlane_distance(r, u, coincident);    break;
    case SurfaceType::SurfaceZPlane:    return SurfaceZPlane_distance(r, u, coincident);    break;
    case SurfaceType::SurfacePlane:     return SurfacePlane_distance(r, u, coincident);     break;
    case SurfaceType::SurfaceXCylinder: return SurfaceXCylinder_distance(r, u, coincident); break;
    case SurfaceType::SurfaceYCylinder: return SurfaceYCylinder_distance(r, u, coincident); break;
    case SurfaceType::SurfaceZCylinder: return SurfaceZCylinder_distance(r, u, coincident); break;
    case SurfaceType::SurfaceSphere:    return SurfaceSphere_distance(r, u, coincident);    break;
    case SurfaceType::SurfaceXCone:     return SurfaceXCone_distance(r, u, coincident);     break;
    case SurfaceType::SurfaceYCone:     return SurfaceYCone_distance(r, u, coincident);     break;
    case SurfaceType::SurfaceZCone:     return SurfaceZCone_distance(r, u, coincident);     break;
    case SurfaceType::SurfaceQuadric:   return SurfaceQuadric_distance(r, u, coincident);   break;
  }
}

Direction Surface::normal(Position r) const
{
  switch(type_){
    case SurfaceType::SurfaceXPlane:    return SurfaceXPlane_normal(r);    break;
    case SurfaceType::SurfaceYPlane:    return SurfaceYPlane_normal(r);    break;
    case SurfaceType::SurfaceZPlane:    return SurfaceZPlane_normal(r);    break;
    case SurfaceType::SurfacePlane:     return SurfacePlane_normal(r);     break;
    case SurfaceType::SurfaceXCylinder: return SurfaceXCylinder_normal(r); break;
    case SurfaceType::SurfaceYCylinder: return SurfaceYCylinder_normal(r); break;
    case SurfaceType::SurfaceZCylinder: return SurfaceZCylinder_normal(r); break;
    case SurfaceType::SurfaceSphere:    return SurfaceSphere_normal(r);    break;
    case SurfaceType::SurfaceXCone:     return SurfaceXCone_normal(r);     break;
    case SurfaceType::SurfaceYCone:     return SurfaceYCone_normal(r);     break;
    case SurfaceType::SurfaceZCone:     return SurfaceZCone_normal(r);     break;
    case SurfaceType::SurfaceQuadric:   return SurfaceQuadric_normal(r);   break;
  }
}

void Surface::to_hdf5_inner(hid_t group_id) const
{
  switch(type_){
    case SurfaceType::SurfaceXPlane:    SurfaceXPlane_to_hdf5_inner(group_id);    break;
    case SurfaceType::SurfaceYPlane:    SurfaceYPlane_to_hdf5_inner(group_id);    break;
    case SurfaceType::SurfaceZPlane:    SurfaceZPlane_to_hdf5_inner(group_id);    break;
    case SurfaceType::SurfacePlane:     SurfacePlane_to_hdf5_inner(group_id);     break;
    case SurfaceType::SurfaceXCylinder: SurfaceXCylinder_to_hdf5_inner(group_id); break;
    case SurfaceType::SurfaceYCylinder: SurfaceYCylinder_to_hdf5_inner(group_id); break;
    case SurfaceType::SurfaceZCylinder: SurfaceZCylinder_to_hdf5_inner(group_id); break;
    case SurfaceType::SurfaceSphere:    SurfaceSphere_to_hdf5_inner(group_id);    break;
    case SurfaceType::SurfaceXCone:     SurfaceXCone_to_hdf5_inner(group_id);     break;
    case SurfaceType::SurfaceYCone:     SurfaceYCone_to_hdf5_inner(group_id);     break;
    case SurfaceType::SurfaceZCone:     SurfaceZCone_to_hdf5_inner(group_id);     break;
    case SurfaceType::SurfaceQuadric:   SurfaceQuadric_to_hdf5_inner(group_id);   break;
  }
}

// Not all types have the ones below...
BoundingBox Surface::bounding_box(bool pos_side) const
{
  switch(type_){
    case SurfaceType::SurfaceXPlane:    return SurfaceXPlane_bounding_box(pos_side);    break;
    case SurfaceType::SurfaceYPlane:    return SurfaceYPlane_bounding_box(pos_side);    break;
    case SurfaceType::SurfaceZPlane:    return SurfaceZPlane_bounding_box(pos_side);    break;
    //case SurfaceType::SurfacePlane:     return SurfacePlane_bounding_box(pos_side);     break;
    case SurfaceType::SurfaceXCylinder: return SurfaceXCylinder_bounding_box(pos_side); break;
    case SurfaceType::SurfaceYCylinder: return SurfaceYCylinder_bounding_box(pos_side); break;
    case SurfaceType::SurfaceZCylinder: return SurfaceZCylinder_bounding_box(pos_side); break;
    case SurfaceType::SurfaceSphere:    return SurfaceSphere_bounding_box(pos_side);    break;
    //case SurfaceType::SurfaceXCone:     return SurfaceXCone_bounding_box(pos_side);     break;
    //case SurfaceType::SurfaceYCone:     return SurfaceYCone_bounding_box(pos_side);     break;
    //case SurfaceType::SurfaceZCone:     return SurfaceZCone_bounding_box(pos_side);     break;
    //case SurfaceType::SurfaceQuadric:   return SurfaceQuadric_bounding_box(pos_side);   break;
    default: return {};
  }
}

bool Surface::periodic_translate(const Surface* other, Position& r,
                                      Direction& u) const
{
  switch(type_){
    case SurfaceType::SurfaceXPlane:    return SurfaceXPlane_periodic_translate(other, r, u);    break;
    case SurfaceType::SurfaceYPlane:    return SurfaceYPlane_periodic_translate(other, r, u);    break;
    case SurfaceType::SurfaceZPlane:    return SurfaceZPlane_periodic_translate(other, r, u);    break;
    case SurfaceType::SurfacePlane:     return SurfacePlane_periodic_translate(other, r, u);     break;
                           /*
    case SurfaceType::SurfaceXCylinder: return SurfaceXCylinder_periodic_translate(other, r, u); break;
    case SurfaceType::SurfaceYCylinder: return SurfaceYCylinder_periodic_translate(other, r, u); break;
    case SurfaceType::SurfaceZCylinder: return SurfaceZCylinder_periodic_translate(other, r, u); break;
    case SurfaceType::SurfaceSphere:    return SurfaceSphere_periodic_translate(other, r, u);    break;
    case SurfaceType::SurfaceXCone:     return SurfaceXCone_periodic_translate(other, r, u);     break;
    case SurfaceType::SurfaceYCone:     return SurfaceYCone_periodic_translate(other, r, u);     break;
    case SurfaceType::SurfaceZCone:     return SurfaceZCone_periodic_translate(other, r, u);     break;
    case SurfaceType::SurfaceQuadric:   return SurfaceQuadric_periodic_translate(other, r, u);   break;
    */
    default: return false;
  }
}

//==============================================================================

void read_surfaces(pugi::xml_node node)
{
  // Count the number of surfaces.
  int n_surfaces = 0;
  for (pugi::xml_node surf_node : node.children("surface")) {n_surfaces++;}
  if (n_surfaces == 0) {
    fatal_error("No surfaces found in geometry.xml!");
  }

  // Loop over XML surface elements and populate the array.
  model::surfaces.reserve(n_surfaces);
  {
    pugi::xml_node surf_node;
    int i_surf;
    for (surf_node = node.child("surface"), i_surf = 0; surf_node;
         surf_node = surf_node.next_sibling("surface"), i_surf++) {
      std::string surf_type = get_node_value(surf_node, "type", true, true);

      /*
      if (surf_type == "x-plane") {
        model::surfaces.push_back(std::make_unique<Surface>(surf_node, Surface::SurfaceType::SurfaceXPlane));

      } else if (surf_type == "y-plane") {
        model::surfaces.push_back(std::make_unique<Surface>(surf_node, Surface::SurfaceType::SurfaceYPlane));

      } else if (surf_type == "z-plane") {
        model::surfaces.push_back(std::make_unique<Surface>(surf_node, Surface::SurfaceType::SurfaceZPlane));

      } else if (surf_type == "plane") {
        model::surfaces.push_back(std::make_unique<Surface>(surf_node, Surface::SurfaceType::SurfacePlane));

      } else if (surf_type == "x-cylinder") {
        model::surfaces.push_back(std::make_unique<Surface>(surf_node, Surface::SurfaceType::SurfaceXCylinder));

      } else if (surf_type == "y-cylinder") {
        model::surfaces.push_back(std::make_unique<Surface>(surf_node, Surface::SurfaceType::SurfaceYCylinder));

      } else if (surf_type == "z-cylinder") {
        model::surfaces.push_back(std::make_unique<Surface>(surf_node, Surface::SurfaceType::SurfaceZCylinder));

      } else if (surf_type == "sphere") {
        model::surfaces.push_back(std::make_unique<Surface>(surf_node, Surface::SurfaceType::SurfaceSphere));

      } else if (surf_type == "x-cone") {
        model::surfaces.push_back(std::make_unique<Surface>(surf_node, Surface::SurfaceType::SurfaceXCone));

      } else if (surf_type == "y-cone") {
        model::surfaces.push_back(std::make_unique<Surface>(surf_node, Surface::SurfaceType::SurfaceYCone));

      } else if (surf_type == "z-cone") {
        model::surfaces.push_back(std::make_unique<Surface>(surf_node, Surface::SurfaceType::SurfaceZCone));

      } else if (surf_type == "quadric") {
        model::surfaces.push_back(std::make_unique<Surface>(surf_node, Surface::SurfaceType::SurfaceQuadric));

      } else {
        fatal_error(fmt::format("Invalid surface type, \"{}\"", surf_type));
      }
      */

      if (surf_type == "x-plane") {
        model::surfaces.emplace_back(surf_node, Surface::SurfaceType::SurfaceXPlane);

      } else if (surf_type == "y-plane") {
        model::surfaces.emplace_back(surf_node, Surface::SurfaceType::SurfaceYPlane);

      } else if (surf_type == "z-plane") {
        model::surfaces.emplace_back(surf_node, Surface::SurfaceType::SurfaceZPlane);

      } else if (surf_type == "plane") {
        model::surfaces.emplace_back(surf_node, Surface::SurfaceType::SurfacePlane);

      } else if (surf_type == "x-cylinder") {
        model::surfaces.emplace_back(surf_node, Surface::SurfaceType::SurfaceXCylinder);

      } else if (surf_type == "y-cylinder") {
        model::surfaces.emplace_back(surf_node, Surface::SurfaceType::SurfaceYCylinder);

      } else if (surf_type == "z-cylinder") {
        model::surfaces.emplace_back(surf_node, Surface::SurfaceType::SurfaceZCylinder);

      } else if (surf_type == "sphere") {
        model::surfaces.emplace_back(surf_node, Surface::SurfaceType::SurfaceSphere);

      } else if (surf_type == "x-cone") {
        model::surfaces.emplace_back(surf_node, Surface::SurfaceType::SurfaceXCone);

      } else if (surf_type == "y-cone") {
        model::surfaces.emplace_back(surf_node, Surface::SurfaceType::SurfaceYCone);

      } else if (surf_type == "z-cone") {
        model::surfaces.emplace_back(surf_node, Surface::SurfaceType::SurfaceZCone);

      } else if (surf_type == "quadric") {
        model::surfaces.emplace_back(surf_node, Surface::SurfaceType::SurfaceQuadric);

      } else {
        fatal_error(fmt::format("Invalid surface type, \"{}\"", surf_type));
      }
    }
  }

  // Fill the surface map.
  for (int i_surf = 0; i_surf < model::surfaces.size(); i_surf++) {
    int id = model::surfaces[i_surf].id_;
    auto in_map = model::surface_map.find(id);
    if (in_map == model::surface_map.end()) {
      model::surface_map[id] = i_surf;
    } else {
      fatal_error(fmt::format(
        "Two or more surfaces use the same unique ID: {}", id));
    }
  }

  // Find the global bounding box (of periodic BC surfaces).
  double xmin {INFTY}, xmax {-INFTY}, ymin {INFTY}, ymax {-INFTY},
         zmin {INFTY}, zmax {-INFTY};
  int i_xmin, i_xmax, i_ymin, i_ymax, i_zmin, i_zmax;
  for (int i_surf = 0; i_surf < model::surfaces.size(); i_surf++) {
    if (model::surfaces[i_surf].bc_ == Surface::BoundaryType::PERIODIC) {
      // Downcast to the PeriodicSurface type.
      Surface* surf_base = &model::surfaces[i_surf];
      //auto surf = dynamic_cast<PeriodicSurface*>(surf_base);
      auto surf = surf_base;

      // Make sure this surface inherits from PeriodicSurface.
      //if (!surf) {
      if (
          surf->type_ != Surface::SurfaceType::SurfaceXPlane &&
          surf->type_ != Surface::SurfaceType::SurfaceYPlane &&
          surf->type_ != Surface::SurfaceType::SurfaceZPlane &&
          surf->type_ != Surface::SurfaceType::SurfacePlane
          ) {
        fatal_error(fmt::format(
          "Periodic boundary condition not supported for surface {}. Periodic "
          "BCs are only supported for planar surfaces.", surf_base->id_));
      }

      // See if this surface makes part of the global bounding box.
      auto bb = surf->bounding_box(true) & surf->bounding_box(false);
      if (bb.xmin > -INFTY && bb.xmin < xmin) {
        xmin = bb.xmin;
        i_xmin = i_surf;
      }
      if (bb.xmax < INFTY && bb.xmax > xmax) {
        xmax = bb.xmax;
        i_xmax = i_surf;
      }
      if (bb.ymin > -INFTY && bb.ymin < ymin) {
        ymin = bb.ymin;
        i_ymin = i_surf;
      }
      if (bb.ymax < INFTY && bb.ymax > ymax) {
        ymax = bb.ymax;
        i_ymax = i_surf;
      }
      if (bb.zmin > -INFTY && bb.zmin < zmin) {
        zmin = bb.zmin;
        i_zmin = i_surf;
      }
      if (bb.zmax < INFTY && bb.zmax > zmax) {
        zmax = bb.zmax;
        i_zmax = i_surf;
      }
    }
  }

  // Set i_periodic for periodic BC surfaces.
  for (int i_surf = 0; i_surf < model::surfaces.size(); i_surf++) {
    if (model::surfaces[i_surf].bc_ == Surface::BoundaryType::PERIODIC) {
      // Downcast to the PeriodicSurface type.
      Surface* surf_base = &model::surfaces[i_surf];
      //auto surf = dynamic_cast<PeriodicSurface*>(surf_base);
      auto surf = surf_base;

      // Also try downcasting to the SurfacePlane type (which must be handled
      // differently).
      //SurfacePlane* surf_p = dynamic_cast<SurfacePlane*>(surf);
      auto surf_p = surf;

      //if (!surf_p) {
      if (surf_p->type_ != Surface::SurfaceType::SurfacePlane) {
        // This is not a SurfacePlane.
        if (surf->i_periodic_ == C_NONE) {
          // The user did not specify the matching periodic surface.  See if we
          // can find the partnered surface from the bounding box information.
          if (i_surf == i_xmin) {
            surf->i_periodic_ = i_xmax;
          } else if (i_surf == i_xmax) {
            surf->i_periodic_ = i_xmin;
          } else if (i_surf == i_ymin) {
            surf->i_periodic_ = i_ymax;
          } else if (i_surf == i_ymax) {
            surf->i_periodic_ = i_ymin;
          } else if (i_surf == i_zmin) {
            surf->i_periodic_ = i_zmax;
          } else if (i_surf == i_zmax) {
            surf->i_periodic_ = i_zmin;
          } else {
            fatal_error("Periodic boundary condition applied to interior "
                        "surface");
          }
        } else {
          // Convert the surface id to an index.
          surf->i_periodic_ = model::surface_map[surf->i_periodic_];
        }
      } else {
        // This is a SurfacePlane.  We won't try to find it's partner if the
        // user didn't specify one.
        if (surf->i_periodic_ == C_NONE) {
          fatal_error(fmt::format("No matching periodic surface specified for "
            "periodic boundary condition on surface {}", surf->id_));
        } else {
          // Convert the surface id to an index.
          surf->i_periodic_ = model::surface_map[surf->i_periodic_];
        }
      }

      // Make sure the opposite surface is also periodic.
      if (model::surfaces[surf->i_periodic_].bc_ != Surface::BoundaryType::PERIODIC) {
        fatal_error(fmt::format("Could not find matching surface for periodic "
          "boundary condition on surface {}", surf->id_));
      }
    }
  }

  // Check to make sure a boundary condition was applied to at least one
  // surface
  bool boundary_exists = false;
  for (const auto& surf : model::surfaces) {
    if (surf.bc_ != Surface::BoundaryType::TRANSMIT) {
      boundary_exists = true;
      break;
    }
  }
  if (settings::run_mode != RunMode::PLOTTING && !boundary_exists) {
    fatal_error("No boundary conditions were applied to any surfaces!");
  }
}

void free_memory_surfaces()
{
  model::surfaces.clear();
  model::surface_map.clear();
}

} // namespace openmc
