#include "openmc/surface.h"

#include <array>
#include <cmath>
#include <sstream>
#include <utility>

#include "openmc/error.h"
#include "openmc/hdf5_interface.h"
#include "openmc/xml_interface.h"


namespace openmc {

//==============================================================================
// Module constant definitions
//==============================================================================

extern "C" const int BC_TRANSMIT {0};
extern "C" const int BC_VACUUM {1};
extern "C" const int BC_REFLECT {2};
extern "C" const int BC_PERIODIC {3};

//==============================================================================
// Global variables
//==============================================================================

int32_t n_surfaces;

std::vector<Surface*> global_surfaces;

std::map<int, int> surface_map;

//==============================================================================
// Helper functions for reading the "coeffs" node of an XML surface element
//==============================================================================

int word_count(const std::string &text)
{
  bool in_word = false;
  int count {0};
  for (auto c = text.begin(); c != text.end(); c++) {
    if (std::isspace(*c)) {
      if (in_word) {
        in_word = false;
        count++;
      }
    } else {
      in_word = true;
    }
  }
  if (in_word) count++;
  return count;
}

void read_coeffs(pugi::xml_node surf_node, int surf_id, double &c1)
{
  // Check the given number of coefficients.
  std::string coeffs = get_node_value(surf_node, "coeffs");
  int n_words = word_count(coeffs);
  if (n_words != 1) {
    std::stringstream err_msg;
    err_msg << "Surface " << surf_id << " expects 1 coeff but was given "
            << n_words;
    fatal_error(err_msg);
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
    std::stringstream err_msg;
    err_msg << "Surface " << surf_id << " expects 3 coeffs but was given "
            << n_words;
    fatal_error(err_msg);
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
    std::stringstream err_msg;
    err_msg << "Surface " << surf_id << " expects 4 coeffs but was given "
            << n_words;
    fatal_error(err_msg);
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
    std::stringstream err_msg;
    err_msg << "Surface " << surf_id << " expects 10 coeffs but was given "
            << n_words;
    fatal_error(err_msg);
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

Surface::Surface(pugi::xml_node surf_node)
{
  if (check_for_node(surf_node, "id")) {
    id = std::stoi(get_node_value(surf_node, "id"));
  } else {
    fatal_error("Must specify id of surface in geometry XML file.");
  }

  if (check_for_node(surf_node, "name")) {
    name = get_node_value(surf_node, "name", false);
  }

  if (check_for_node(surf_node, "boundary")) {
    std::string surf_bc = get_node_value(surf_node, "boundary", true, true);

    if (surf_bc == "transmission" || surf_bc == "transmit" ||surf_bc.empty()) {
      bc = BC_TRANSMIT;

    } else if (surf_bc == "vacuum") {
      bc = BC_VACUUM;

    } else if (surf_bc == "reflective" || surf_bc == "reflect"
               || surf_bc == "reflecting") {
      bc = BC_REFLECT;
    } else if (surf_bc == "periodic") {
      bc = BC_PERIODIC;
    } else {
      std::stringstream err_msg;
      err_msg << "Unknown boundary condition \"" << surf_bc
              << "\" specified on surface " << id;
      fatal_error(err_msg);
    }

  } else {
    bc = BC_TRANSMIT;
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
Surface::reflect(Position r, Direction u) const
{
  // Determine projection of direction onto normal and squared magnitude of
  // normal.
  Direction n = normal(r);
  const double projection = n.dot(u);
  const double magnitude = n.dot(n);

  // Reflect direction according to normal.
  return u -= (2.0 * projection / magnitude) * n;
}

void
Surface::to_hdf5(hid_t group_id) const
{
  std::string group_name {"surface "};
  group_name += std::to_string(id);

  hid_t surf_group = create_group(group_id, group_name);

  switch(bc) {
    case BC_TRANSMIT :
      write_string(surf_group, "boundary_type", "transmission", false);
      break;
    case BC_VACUUM :
      write_string(surf_group, "boundary_type", "vacuum", false);
      break;
    case BC_REFLECT :
      write_string(surf_group, "boundary_type", "reflective", false);
      break;
    case BC_PERIODIC :
      write_string(surf_group, "boundary_type", "periodic", false);
      break;
  }

  if (!name.empty()) {
    write_string(surf_group, "name", name, false);
  }

  to_hdf5_inner(surf_group);

  close_group(surf_group);
}

//==============================================================================
// PeriodicSurface implementation
//==============================================================================

PeriodicSurface::PeriodicSurface(pugi::xml_node surf_node)
  : Surface {surf_node}
{
  if (check_for_node(surf_node, "periodic_surface_id")) {
    i_periodic = std::stoi(get_node_value(surf_node, "periodic_surface_id"));
  }
}

//==============================================================================
// Generic functions for x-, y-, and z-, planes.
//==============================================================================

// The template parameter indicates the axis normal to the plane.
template<int i> double
axis_aligned_plane_distance(Position r, Direction u, bool coincident, double offset)
{
  const double f = offset - r[i];
  if (coincident or std::abs(f) < FP_COINCIDENT or u[i] == 0.0) return INFTY;
  const double d = f / u[i];
  if (d < 0.0) return INFTY;
  return d;
}

//==============================================================================
// SurfaceXPlane implementation
//==============================================================================

SurfaceXPlane::SurfaceXPlane(pugi::xml_node surf_node)
  : PeriodicSurface(surf_node)
{
  read_coeffs(surf_node, id, x0);
}

double SurfaceXPlane::evaluate(Position r) const
{
  return r.x - x0;
}

double SurfaceXPlane::distance(Position r, Direction u, bool coincident) const
{
  return axis_aligned_plane_distance<0>(r, u, coincident, x0);
}

Direction SurfaceXPlane::normal(Position r) const
{
  return {1., 0., 0.};
}

void SurfaceXPlane::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "x-plane", false);
  std::array<double, 1> coeffs {{x0}};
  write_dataset(group_id, "coefficients", coeffs);
}

bool SurfaceXPlane::periodic_translate(const PeriodicSurface* other,
                                       Position& r,  Direction& u) const
{
  Direction other_n = other->normal(r);
  if (other_n.x == 1 and other_n.y == 0 and other_n.z == 0) {
    r.x = x0;
    return false;
  } else {
    // Assume the partner is an YPlane (the only supported partner).  Use the
    // evaluate function to find y0, then adjust position/Direction for rotational
    // symmetry.
    double y0 = -other->evaluate({0., 0., 0.});
    r.y = r.x - x0 + y0;
    r.x = x0;

    double ux = u.x;
    u.x = -u.y;
    u.y = ux;

    return true;
  }
}

BoundingBox
SurfaceXPlane::bounding_box() const
{
  return {x0, x0, -INFTY, INFTY, -INFTY, INFTY};
}

//==============================================================================
// SurfaceYPlane implementation
//==============================================================================

SurfaceYPlane::SurfaceYPlane(pugi::xml_node surf_node)
  : PeriodicSurface(surf_node)
{
  read_coeffs(surf_node, id, y0);
}

double SurfaceYPlane::evaluate(Position r) const
{
  return r.y - y0;
}

double SurfaceYPlane::distance(Position r, Direction u, bool coincident) const
{
  return axis_aligned_plane_distance<1>(r, u, coincident, y0);
}

Direction SurfaceYPlane::normal(Position r) const
{
  return {0., 1., 0.};
}

void SurfaceYPlane::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "y-plane", false);
  std::array<double, 1> coeffs {{y0}};
  write_dataset(group_id, "coefficients", coeffs);
}

bool SurfaceYPlane::periodic_translate(const PeriodicSurface *other,
                                       Position& r, Direction& u) const
{
  Direction other_n = other->normal(r);
  if (other_n.x == 0 and other_n.y == 1 and other_n.z == 0) {
    // The periodic partner is also aligned along y.  Just change the y coord.
    r.y = y0;
    return false;
  } else {
    // Assume the partner is an XPlane (the only supported partner).  Use the
    // evaluate function to find x0, then adjust position/Direction for rotational
    // symmetry.
    double x0 = -other->evaluate({0., 0., 0.});
    r.x = r.y - y0 + x0;
    r.y = y0;

    double ux = u.x;
    u.x = u.y;
    u.y = -ux;

    return true;
  }
}

BoundingBox
SurfaceYPlane::bounding_box() const
{
  return {-INFTY, INFTY, y0, y0, -INFTY, INFTY};
}

//==============================================================================
// SurfaceZPlane implementation
//==============================================================================

SurfaceZPlane::SurfaceZPlane(pugi::xml_node surf_node)
  : PeriodicSurface(surf_node)
{
  read_coeffs(surf_node, id, z0);
}

double SurfaceZPlane::evaluate(Position r) const
{
  return r.z - z0;
}

double SurfaceZPlane::distance(Position r, Direction u, bool coincident) const
{
  return axis_aligned_plane_distance<2>(r, u, coincident, z0);
}

Direction SurfaceZPlane::normal(Position r) const
{
  return {0., 0., 1.};
}

void SurfaceZPlane::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "z-plane", false);
  std::array<double, 1> coeffs {{z0}};
  write_dataset(group_id, "coefficients", coeffs);
}

bool SurfaceZPlane::periodic_translate(const PeriodicSurface* other,
                                       Position& r, Direction& u) const
{
  // Assume the other plane is aligned along z.  Just change the z coord.
  r.z = z0;
  return false;
}

BoundingBox
SurfaceZPlane::bounding_box() const
{
  return {-INFTY, INFTY, -INFTY, INFTY, z0, z0};
}

//==============================================================================
// SurfacePlane implementation
//==============================================================================

SurfacePlane::SurfacePlane(pugi::xml_node surf_node)
  : PeriodicSurface(surf_node)
{
  read_coeffs(surf_node, id, A, B, C, D);
}

double
SurfacePlane::evaluate(Position r) const
{
  return A*r.x + B*r.y + C*r.z - D;
}

double
SurfacePlane::distance(Position r, Direction u, bool coincident) const
{
  const double f = A*r.x + B*r.y + C*r.z - D;
  const double projection = A*u.x + B*u.y + C*u.z;
  if (coincident or std::abs(f) < FP_COINCIDENT or projection == 0.0) {
    return INFTY;
  } else {
    const double d = -f / projection;
    if (d < 0.0) return INFTY;
    return d;
  }
}

Direction
SurfacePlane::normal(Position r) const
{
  return {A, B, C};
}

void SurfacePlane::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "plane", false);
  std::array<double, 4> coeffs {{A, B, C, D}};
  write_dataset(group_id, "coefficients", coeffs);
}

bool SurfacePlane::periodic_translate(const PeriodicSurface* other, Position& r,
                                      Direction& u) const
{
  // This function assumes the other plane shares this plane's normal direction.

  // Determine the distance to intersection.
  double d = evaluate(r) / (A*A + B*B + C*C);

  // Move the particle that distance along the normal vector.
  r.x -= d * A;
  r.y -= d * B;
  r.z -= d * C;

  return false;
}

BoundingBox
SurfacePlane::bounding_box() const
{
  return {-INFTY, INFTY, -INFTY, INFTY, -INFTY, INFTY};
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

  } else if (coincident or std::abs(c) < FP_COINCIDENT) {
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

SurfaceXCylinder::SurfaceXCylinder(pugi::xml_node surf_node)
  : Surface(surf_node)
{
  read_coeffs(surf_node, id, y0, z0, radius);
}

double SurfaceXCylinder::evaluate(Position r) const
{
  return axis_aligned_cylinder_evaluate<1, 2>(r, y0, z0, radius);
}

double SurfaceXCylinder::distance(Position r, Direction u, bool coincident) const
{
  return axis_aligned_cylinder_distance<0, 1, 2>(r, u, coincident, y0, z0,
                                                 radius);
}

Direction SurfaceXCylinder::normal(Position r) const
{
  return axis_aligned_cylinder_normal<0, 1, 2>(r, y0, z0);
}


void SurfaceXCylinder::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "x-cylinder", false);
  std::array<double, 3> coeffs {{y0, z0, radius}};
  write_dataset(group_id, "coefficients", coeffs);
}

//==============================================================================
// SurfaceYCylinder implementation
//==============================================================================

SurfaceYCylinder::SurfaceYCylinder(pugi::xml_node surf_node)
  : Surface(surf_node)
{
  read_coeffs(surf_node, id, x0, z0, radius);
}

double SurfaceYCylinder::evaluate(Position r) const
{
  return axis_aligned_cylinder_evaluate<0, 2>(r, x0, z0, radius);
}

double SurfaceYCylinder::distance(Position r, Direction u, bool coincident) const
{
  return axis_aligned_cylinder_distance<1, 0, 2>(r, u, coincident, x0, z0,
                                                 radius);
}

Direction SurfaceYCylinder::normal(Position r) const
{
  return axis_aligned_cylinder_normal<1, 0, 2>(r, x0, z0);
}

void SurfaceYCylinder::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "y-cylinder", false);
  std::array<double, 3> coeffs {{x0, z0, radius}};
  write_dataset(group_id, "coefficients", coeffs);
}

//==============================================================================
// SurfaceZCylinder implementation
//==============================================================================

SurfaceZCylinder::SurfaceZCylinder(pugi::xml_node surf_node)
  : Surface(surf_node)
{
  read_coeffs(surf_node, id, x0, y0, radius);
}

double SurfaceZCylinder::evaluate(Position r) const
{
  return axis_aligned_cylinder_evaluate<0, 1>(r, x0, y0, radius);
}

double SurfaceZCylinder::distance(Position r, Direction u, bool coincident) const
{
  return axis_aligned_cylinder_distance<2, 0, 1>(r, u, coincident, x0, y0,
                                                 radius);
}

Direction SurfaceZCylinder::normal(Position r) const
{
  return axis_aligned_cylinder_normal<2, 0, 1>(r, x0, y0);
}

void SurfaceZCylinder::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "z-cylinder", false);
  std::array<double, 3> coeffs {{x0, y0, radius}};
  write_dataset(group_id, "coefficients", coeffs);
}

//==============================================================================
// SurfaceSphere implementation
//==============================================================================

SurfaceSphere::SurfaceSphere(pugi::xml_node surf_node)
  : Surface(surf_node)
{
  read_coeffs(surf_node, id, x0, y0, z0, radius);
}

double SurfaceSphere::evaluate(Position r) const
{
  const double x = r.x - x0;
  const double y = r.y - y0;
  const double z = r.z - z0;
  return x*x + y*y + z*z - radius*radius;
}

double SurfaceSphere::distance(Position r, Direction u, bool coincident) const
{
  const double x = r.x - x0;
  const double y = r.y - y0;
  const double z = r.z - z0;
  const double k = x*u.x + y*u.y + z*u.z;
  const double c = x*x + y*y + z*z - radius*radius;
  const double quad = k*k - c;

  if (quad < 0.0) {
    // No intersection with sphere.
    return INFTY;

  } else if (coincident or std::abs(c) < FP_COINCIDENT) {
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

Direction SurfaceSphere::normal(Position r) const
{
  return {2.0*(r.x - x0), 2.0*(r.y - y0), 2.0*(r.z - z0)};
}

void SurfaceSphere::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "sphere", false);
  std::array<double, 4> coeffs {{x0, y0, z0, radius}};
  write_dataset(group_id, "coefficients", coeffs);
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

  } else if (coincident or std::abs(c) < FP_COINCIDENT) {
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

SurfaceXCone::SurfaceXCone(pugi::xml_node surf_node)
  : Surface(surf_node)
{
  read_coeffs(surf_node, id, x0, y0, z0, radius_sq);
}

double SurfaceXCone::evaluate(Position r) const
{
  return axis_aligned_cone_evaluate<0, 1, 2>(r, x0, y0, z0, radius_sq);
}

double SurfaceXCone::distance(Position r, Direction u, bool coincident) const
{
  return axis_aligned_cone_distance<0, 1, 2>(r, u, coincident, x0, y0, z0,
                                             radius_sq);
}

Direction SurfaceXCone::normal(Position r) const
{
  return axis_aligned_cone_normal<0, 1, 2>(r, x0, y0, z0, radius_sq);
}

void SurfaceXCone::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "x-cone", false);
  std::array<double, 4> coeffs {{x0, y0, z0, radius_sq}};
  write_dataset(group_id, "coefficients", coeffs);
}

//==============================================================================
// SurfaceYCone implementation
//==============================================================================

SurfaceYCone::SurfaceYCone(pugi::xml_node surf_node)
  : Surface(surf_node)
{
  read_coeffs(surf_node, id, x0, y0, z0, radius_sq);
}

double SurfaceYCone::evaluate(Position r) const
{
  return axis_aligned_cone_evaluate<1, 0, 2>(r, y0, x0, z0, radius_sq);
}

double SurfaceYCone::distance(Position r, Direction u, bool coincident) const
{
  return axis_aligned_cone_distance<1, 0, 2>(r, u, coincident, y0, x0, z0,
                                             radius_sq);
}

Direction SurfaceYCone::normal(Position r) const
{
  return axis_aligned_cone_normal<1, 0, 2>(r, y0, x0, z0, radius_sq);
}

void SurfaceYCone::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "y-cone", false);
  std::array<double, 4> coeffs {{x0, y0, z0, radius_sq}};
  write_dataset(group_id, "coefficients", coeffs);
}

//==============================================================================
// SurfaceZCone implementation
//==============================================================================

SurfaceZCone::SurfaceZCone(pugi::xml_node surf_node)
  : Surface(surf_node)
{
  read_coeffs(surf_node, id, x0, y0, z0, radius_sq);
}

double SurfaceZCone::evaluate(Position r) const
{
  return axis_aligned_cone_evaluate<2, 0, 1>(r, z0, x0, y0, radius_sq);
}

double SurfaceZCone::distance(Position r, Direction u, bool coincident) const
{
  return axis_aligned_cone_distance<2, 0, 1>(r, u, coincident, z0, x0, y0,
                                             radius_sq);
}

Direction SurfaceZCone::normal(Position r) const
{
  return axis_aligned_cone_normal<2, 0, 1>(r, z0, x0, y0, radius_sq);
}

void SurfaceZCone::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "z-cone", false);
  std::array<double, 4> coeffs {{x0, y0, z0, radius_sq}};
  write_dataset(group_id, "coefficients", coeffs);
}

//==============================================================================
// SurfaceQuadric implementation
//==============================================================================

SurfaceQuadric::SurfaceQuadric(pugi::xml_node surf_node)
  : Surface(surf_node)
{
  read_coeffs(surf_node, id, A, B, C, D, E, F, G, H, J, K);
}

double
SurfaceQuadric::evaluate(Position r) const
{
  const double x = r.x;
  const double y = r.y;
  const double z = r.z;
  return x*(A*x + D*y + G) +
         y*(B*y + E*z + H) +
         z*(C*z + F*x + J) + K;
}

double
SurfaceQuadric::distance(Position r, Direction ang, bool coincident) const
{
  const double &x = r.x;
  const double &y = r.y;
  const double &z = r.z;
  const double &u = ang.x;
  const double &v = ang.y;
  const double &w = ang.z;

  const double a = A*u*u + B*v*v + C*w*w + D*u*v + E*v*w + F*u*w;
  const double k = (A*u*x + B*v*y + C*w*z + 0.5*(D*(u*y + v*x) +
                    E*(v*z + w*y) + F*(w*x + u*z) + G*u + H*v + J*w));
  const double c = A*x*x + B*y*y + C*z*z + D*x*y + E*y*z +  F*x*z + G*x + H*y
                   + J*z + K;
  double quad = k*k - a*c;

  double d;

  if (quad < 0.0) {
    // No intersection with surface.
    return INFTY;

  } else if (coincident or std::abs(c) < FP_COINCIDENT) {
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
SurfaceQuadric::normal(Position r) const
{
  const double &x = r.x;
  const double &y = r.y;
  const double &z = r.z;
  return {2.0*A*x + D*y + F*z + G,
      2.0*B*y + D*x + E*z + H,
      2.0*C*z + E*y + F*x + J};
}

void SurfaceQuadric::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "quadric", false);
  std::array<double, 10> coeffs {{A, B, C, D, E, F, G, H, J, K}};
  write_dataset(group_id, "coefficients", coeffs);
}

//==============================================================================

extern "C" void
read_surfaces(pugi::xml_node* node)
{
  // Count the number of surfaces.
  for (pugi::xml_node surf_node: node->children("surface")) {n_surfaces++;}
  if (n_surfaces == 0) {
    fatal_error("No surfaces found in geometry.xml!");
  }

  // Loop over XML surface elements and populate the array.
  global_surfaces.reserve(n_surfaces);
  {
    pugi::xml_node surf_node;
    int i_surf;
    for (surf_node = node->child("surface"), i_surf = 0; surf_node;
         surf_node = surf_node.next_sibling("surface"), i_surf++) {
      std::string surf_type = get_node_value(surf_node, "type", true, true);

      if (surf_type == "x-plane") {
        global_surfaces.push_back(new SurfaceXPlane(surf_node));

      } else if (surf_type == "y-plane") {
        global_surfaces.push_back(new SurfaceYPlane(surf_node));

      } else if (surf_type == "z-plane") {
        global_surfaces.push_back(new SurfaceZPlane(surf_node));

      } else if (surf_type == "plane") {
        global_surfaces.push_back(new SurfacePlane(surf_node));

      } else if (surf_type == "x-cylinder") {
        global_surfaces.push_back(new SurfaceXCylinder(surf_node));

      } else if (surf_type == "y-cylinder") {
        global_surfaces.push_back(new SurfaceYCylinder(surf_node));

      } else if (surf_type == "z-cylinder") {
        global_surfaces.push_back(new SurfaceZCylinder(surf_node));

      } else if (surf_type == "sphere") {
        global_surfaces.push_back(new SurfaceSphere(surf_node));

      } else if (surf_type == "x-cone") {
        global_surfaces.push_back(new SurfaceXCone(surf_node));

      } else if (surf_type == "y-cone") {
        global_surfaces.push_back(new SurfaceYCone(surf_node));

      } else if (surf_type == "z-cone") {
        global_surfaces.push_back(new SurfaceZCone(surf_node));

      } else if (surf_type == "quadric") {
        global_surfaces.push_back(new SurfaceQuadric(surf_node));

      } else {
        std::stringstream err_msg;
        err_msg << "Invalid surface type, \"" << surf_type << "\"";
        fatal_error(err_msg);
      }
    }
  }

  // Fill the surface map.
  for (int i_surf = 0; i_surf < n_surfaces; i_surf++) {
    int id = global_surfaces[i_surf]->id;
    auto in_map = surface_map.find(id);
    if (in_map == surface_map.end()) {
      surface_map[id] = i_surf;
    } else {
      std::stringstream err_msg;
      err_msg << "Two or more surfaces use the same unique ID: " << id;
      fatal_error(err_msg);
    }
  }

  // Find the global bounding box (of periodic BC surfaces).
  double xmin {INFTY}, xmax {-INFTY}, ymin {INFTY}, ymax {-INFTY},
         zmin {INFTY}, zmax {-INFTY};
  int i_xmin, i_xmax, i_ymin, i_ymax, i_zmin, i_zmax;
  for (int i_surf = 0; i_surf < n_surfaces; i_surf++) {
    if (global_surfaces[i_surf]->bc == BC_PERIODIC) {
      // Downcast to the PeriodicSurface type.
      Surface* surf_base = global_surfaces[i_surf];
      PeriodicSurface* surf = dynamic_cast<PeriodicSurface*>(surf_base);

      // Make sure this surface inherits from PeriodicSurface.
      if (!surf) {
        std::stringstream err_msg;
        err_msg << "Periodic boundary condition not supported for surface "
                << surf_base->id
                << ". Periodic BCs are only supported for planar surfaces.";
        fatal_error(err_msg);
      }

      // See if this surface makes part of the global bounding box.
      BoundingBox bb = surf->bounding_box();
      if (bb.xmin > -INFTY and bb.xmin < xmin) {
        xmin = bb.xmin;
        i_xmin = i_surf;
      }
      if (bb.xmax < INFTY and bb.xmax > xmax) {
        xmax = bb.xmax;
        i_xmax = i_surf;
      }
      if (bb.ymin > -INFTY and bb.ymin < ymin) {
        ymin = bb.ymin;
        i_ymin = i_surf;
      }
      if (bb.ymax < INFTY and bb.ymax > ymax) {
        ymax = bb.ymax;
        i_ymax = i_surf;
      }
      if (bb.zmin > -INFTY and bb.zmin < zmin) {
        zmin = bb.zmin;
        i_zmin = i_surf;
      }
      if (bb.zmax < INFTY and bb.zmax > zmax) {
        zmax = bb.zmax;
        i_zmax = i_surf;
      }
    }
  }

  // Set i_periodic for periodic BC surfaces.
  for (int i_surf = 0; i_surf < n_surfaces; i_surf++) {
    if (global_surfaces[i_surf]->bc == BC_PERIODIC) {
      // Downcast to the PeriodicSurface type.
      Surface* surf_base = global_surfaces[i_surf];
      PeriodicSurface* surf = dynamic_cast<PeriodicSurface*>(surf_base);

      // Also try downcasting to the SurfacePlane type (which must be handled
      // differently).
      SurfacePlane* surf_p = dynamic_cast<SurfacePlane*>(surf);

      if (!surf_p) {
        // This is not a SurfacePlane.
        if (surf->i_periodic == C_NONE) {
          // The user did not specify the matching periodic surface.  See if we
          // can find the partnered surface from the bounding box information.
          if (i_surf == i_xmin) {
            surf->i_periodic = i_xmax;
          } else if (i_surf == i_xmax) {
            surf->i_periodic = i_xmin;
          } else if (i_surf == i_ymin) {
            surf->i_periodic = i_ymax;
          } else if (i_surf == i_ymax) {
            surf->i_periodic = i_ymin;
          } else if (i_surf == i_zmin) {
            surf->i_periodic = i_zmax;
          } else if (i_surf == i_zmax) {
            surf->i_periodic = i_zmin;
          } else {
            fatal_error("Periodic boundary condition applied to interior "
                        "surface");
          }
        } else {
          // Convert the surface id to an index.
          surf->i_periodic = surface_map[surf->i_periodic];
        }
      } else {
        // This is a SurfacePlane.  We won't try to find it's partner if the
        // user didn't specify one.
        if (surf->i_periodic == C_NONE) {
          std::stringstream err_msg;
          err_msg << "No matching periodic surface specified for periodic "
                     "boundary condition on surface " << surf->id;
          fatal_error(err_msg);
        } else {
          // Convert the surface id to an index.
          surf->i_periodic = surface_map[surf->i_periodic];
        }
      }

      // Make sure the opposite surface is also periodic.
      if (global_surfaces[surf->i_periodic]->bc != BC_PERIODIC) {
        std::stringstream err_msg;
        err_msg << "Could not find matching surface for periodic boundary "
                   "condition on surface " << surf->id;
        fatal_error(err_msg);
      }
    }
  }
}

//==============================================================================
// Fortran compatibility functions
//==============================================================================

extern "C" {
  Surface* surface_pointer(int surf_ind) {return global_surfaces[surf_ind];}

  int surface_id(Surface* surf) {return surf->id;}

  int surface_bc(Surface* surf) {return surf->bc;}

  void surface_reflect(Surface* surf, double xyz[3], double uvw[3])
  {
    Position r {xyz};
    Direction u {uvw};
    u = surf->reflect(r, u);

    uvw[0] = u.x;
    uvw[1] = u.y;
    uvw[2] = u.z;
  }

  void surface_normal(Surface* surf, double xyz[3], double uvw[3])
  {
    Position r {xyz};
    Direction u = surf->normal(r);
    uvw[0] = u.x;
    uvw[1] = u.y;
    uvw[2] = u.z;
  }

  void surface_to_hdf5(Surface* surf, hid_t group) {surf->to_hdf5(group);}

  int surface_i_periodic(PeriodicSurface* surf) {return surf->i_periodic;}

  bool
  surface_periodic(PeriodicSurface* surf, PeriodicSurface* other, double xyz[3],
                   double uvw[3])
  {
    Position r {xyz};
    Direction u {uvw};
    bool rotational = surf->periodic_translate(other, r, u);

    // Copy back to arrays
    xyz[0] = r.x;
    xyz[1] = r.y;
    xyz[2] = r.z;
    uvw[0] = u.x;
    uvw[1] = u.y;
    uvw[2] = u.z;

    return rotational;
  }

  void free_memory_surfaces_c()
  {
    for (Surface* surf : global_surfaces) {delete surf;}
    global_surfaces.clear();
    n_surfaces = 0;
    surface_map.clear();
  }
}

} // namespace openmc
