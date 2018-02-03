#include "surface.h"

#include <array>
#include <cmath>
#include <sstream>

#include "error.h"
#include "hdf5_interface.h"
#include "xml_interface.h"


namespace openmc {

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
    id = stoi(get_node_value(surf_node, "id"));
  } else {
    fatal_error("Must specify id of surface in geometry XML file.");
  }

  if (check_for_node(surf_node, "name")) {
    name = get_node_value(surf_node, "name");
  }

  if (check_for_node(surf_node, "boundary")) {
    std::string surf_bc = get_node_value(surf_node, "boundary");

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
Surface::sense(const double xyz[3], const double uvw[3]) const
{
  // Evaluate the surface equation at the particle's coordinates to determine
  // which side the particle is on.
  const double f = evaluate(xyz);

  // Check which side of surface the point is on.
  if (std::abs(f) < FP_COINCIDENT) {
    // Particle may be coincident with this surface. To determine the sense, we
    // look at the direction of the particle relative to the surface normal (by
    // default in the positive direction) via their dot product.
    double norm[3];
    normal(xyz, norm);
    return uvw[0] * norm[0] + uvw[1] * norm[1] + uvw[2] * norm[2] > 0.0;
  }
  return f > 0.0;
}

void
Surface::reflect(const double xyz[3], double uvw[3]) const
{
  // Determine projection of direction onto normal and squared magnitude of
  // normal.
  double norm[3];
  normal(xyz, norm);
  const double projection = norm[0]*uvw[0] + norm[1]*uvw[1] + norm[2]*uvw[2];
  const double magnitude = norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2];

  // Reflect direction according to normal.
  uvw[0] -= 2.0 * projection / magnitude * norm[0];
  uvw[1] -= 2.0 * projection / magnitude * norm[1];
  uvw[2] -= 2.0 * projection / magnitude * norm[2];
}

void
Surface::to_hdf5(hid_t group_id) const
{
  std::string group_name {"surface "};
  group_name += std::to_string(id);

  hid_t surf_group = create_group(group_id, group_name);

  switch(bc) {
    case BC_TRANSMIT :
      write_string(surf_group, "boundary_type", "transmission");
      break;
    case BC_VACUUM :
      write_string(surf_group, "boundary_type", "vacuum");
      break;
    case BC_REFLECT :
      write_string(surf_group, "boundary_type", "reflective");
      break;
    case BC_PERIODIC :
      write_string(surf_group, "boundary_type", "periodic");
      break;
  }

  if (!name.empty()) {
    write_string(surf_group, "name", name);
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
    i_periodic = stoi(get_node_value(surf_node, "periodic_surface_id"));
  }
}

//==============================================================================
// Generic functions for x-, y-, and z-, planes.
//==============================================================================

// The template parameter indicates the axis normal to the plane.
template<int i> double
axis_aligned_plane_evaluate(const double xyz[3], double offset)
{
  return xyz[i] - offset;
}

// The template parameter indicates the axis normal to the plane.
template<int i> double
axis_aligned_plane_distance(const double xyz[3], const double uvw[3],
                            bool coincident, double offset)
{
  const double f = offset - xyz[i];
  if (coincident or std::abs(f) < FP_COINCIDENT or uvw[i] == 0.0) return INFTY;
  const double d = f / uvw[i];
  if (d < 0.0) return INFTY;
  return d;
}

// The first template parameter indicates the axis normal to the plane.  The
// other two parameters indicate the other two axes.
template<int i1, int i2, int i3> void
axis_aligned_plane_normal(const double xyz[3], double uvw[3])
{
  uvw[i1] = 1.0;
  uvw[i2] = 0.0;
  uvw[i3] = 0.0;
}

//==============================================================================
// SurfaceXPlane implementation
//==============================================================================

SurfaceXPlane::SurfaceXPlane(pugi::xml_node surf_node)
  : PeriodicSurface(surf_node)
{
  read_coeffs(surf_node, id, x0);
}

inline double SurfaceXPlane::evaluate(const double xyz[3]) const
{
  return axis_aligned_plane_evaluate<0>(xyz, x0);
}

inline double SurfaceXPlane::distance(const double xyz[3], const double uvw[3],
                                      bool coincident) const
{
  return axis_aligned_plane_distance<0>(xyz, uvw, coincident, x0);
}

inline void SurfaceXPlane::normal(const double xyz[3], double uvw[3]) const
{
  axis_aligned_plane_normal<0, 1, 2>(xyz, uvw);
}

void SurfaceXPlane::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "x-plane");
  std::array<double, 1> coeffs {{x0}};
  write_double_1D(group_id, "coefficients", coeffs);
}

bool SurfaceXPlane::periodic_translate(PeriodicSurface *other, double xyz[3],
     double uvw[3]) const
{
  double other_norm[3];
  other->normal(xyz, other_norm);
  if (other_norm[0] == 1 and other_norm[1] == 0 and other_norm[2] == 0) {
    xyz[0] = x0;
    return false;
  } else {
    // Assume the partner is an YPlane (the only supported partner).  Use the
    // evaluate function to find y0, then adjust xyz and uvw for rotational
    // symmetry.
    double xyz_test[3] {0, 0, 0};
    double y0 = -other->evaluate(xyz_test);
    xyz[1] = xyz[0] - x0 + y0;
    xyz[0] = x0;

    double u = uvw[0];
    uvw[0] = -uvw[1];
    uvw[1] = u;

    return true;
  }
}

BoundingBox
SurfaceXPlane::bounding_box() const
{
  BoundingBox out {x0, x0, -INFTY, INFTY, -INFTY, INFTY};
  return out;
}

//==============================================================================
// SurfaceYPlane implementation
//==============================================================================

SurfaceYPlane::SurfaceYPlane(pugi::xml_node surf_node)
  : PeriodicSurface(surf_node)
{
  read_coeffs(surf_node, id, y0);
}

inline double SurfaceYPlane::evaluate(const double xyz[3]) const
{
  return axis_aligned_plane_evaluate<1>(xyz, y0);
}

inline double SurfaceYPlane::distance(const double xyz[3], const double uvw[3],
                                      bool coincident) const
{
  return axis_aligned_plane_distance<1>(xyz, uvw, coincident, y0);
}

inline void SurfaceYPlane::normal(const double xyz[3], double uvw[3]) const
{
  axis_aligned_plane_normal<1, 0, 2>(xyz, uvw);
}

void SurfaceYPlane::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "y-plane");
  std::array<double, 1> coeffs {{y0}};
  write_double_1D(group_id, "coefficients", coeffs);
}

bool SurfaceYPlane::periodic_translate(PeriodicSurface *other, double xyz[3],
                                       double uvw[3]) const
{
  double other_norm[3];
  other->normal(xyz, other_norm);
  if (other_norm[0] == 0 and other_norm[1] == 1 and other_norm[2] == 0) {
    // The periodic partner is also aligned along y.  Just change the y coord.
    xyz[1] = y0;
    return false;
  } else {
    // Assume the partner is an XPlane (the only supported partner).  Use the
    // evaluate function to find x0, then adjust xyz and uvw for rotational
    // symmetry.
    double xyz_test[3] {0, 0, 0};
    double x0 = -other->evaluate(xyz_test);
    xyz[0] = xyz[1] - y0 + x0;
    xyz[1] = y0;

    double u = uvw[0];
    uvw[0] = uvw[1];
    uvw[1] = -u;

    return true;
  }
}

BoundingBox
SurfaceYPlane::bounding_box() const
{
  BoundingBox out {-INFTY, INFTY, y0, y0, -INFTY, INFTY};
  return out;
}

//==============================================================================
// SurfaceZPlane implementation
//==============================================================================

SurfaceZPlane::SurfaceZPlane(pugi::xml_node surf_node)
  : PeriodicSurface(surf_node)
{
  read_coeffs(surf_node, id, z0);
}

inline double SurfaceZPlane::evaluate(const double xyz[3]) const
{
  return axis_aligned_plane_evaluate<2>(xyz, z0);
}

inline double SurfaceZPlane::distance(const double xyz[3], const double uvw[3],
                                      bool coincident) const
{
  return axis_aligned_plane_distance<2>(xyz, uvw, coincident, z0);
}

inline void SurfaceZPlane::normal(const double xyz[3], double uvw[3]) const
{
  axis_aligned_plane_normal<2, 0, 1>(xyz, uvw);
}

void SurfaceZPlane::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "z-plane");
  std::array<double, 1> coeffs {{z0}};
  write_double_1D(group_id, "coefficients", coeffs);
}

bool SurfaceZPlane::periodic_translate(PeriodicSurface *other, double xyz[3],
                                       double uvw[3]) const
{
  // Assume the other plane is aligned along z.  Just change the z coord.
  xyz[2] = z0;
  return false;
}

BoundingBox
SurfaceZPlane::bounding_box() const
{
  BoundingBox out {-INFTY, INFTY, -INFTY, INFTY, z0, z0};
  return out;
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
SurfacePlane::evaluate(const double xyz[3]) const
{
  return A*xyz[0] + B*xyz[1] + C*xyz[2] - D;
}

double
SurfacePlane::distance(const double xyz[3], const double uvw[3],
                       bool coincident) const
{
  const double f = A*xyz[0] + B*xyz[1] + C*xyz[2] - D;
  const double projection = A*uvw[0] + B*uvw[1] + C*uvw[2];
  if (coincident or std::abs(f) < FP_COINCIDENT or projection == 0.0) {
    return INFTY;
  } else {
    const double d = -f / projection;
    if (d < 0.0) return INFTY;
    return d;
  }
}

void
SurfacePlane::normal(const double xyz[3], double uvw[3]) const
{
  uvw[0] = A;
  uvw[1] = B;
  uvw[2] = C;
}

void SurfacePlane::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "plane");
  std::array<double, 4> coeffs {{A, B, C, D}};
  write_double_1D(group_id, "coefficients", coeffs);
}

bool SurfacePlane::periodic_translate(PeriodicSurface *other, double xyz[3],
                                      double uvw[3]) const
{
  // This function assumes the other plane shares this plane's normal direction.

  // Determine the distance to intersection.
  double d = evaluate(xyz) / (A*A + B*B + C*C);

  // Move the particle that distance along the normal vector.
  xyz[0] -= d * A;
  xyz[1] -= d * B;
  xyz[2] -= d * C;

  return false;
}

BoundingBox
SurfacePlane::bounding_box() const
{
  BoundingBox out {-INFTY, INFTY, -INFTY, INFTY, -INFTY, INFTY};
  return out;
}

//==============================================================================
// Generic functions for x-, y-, and z-, cylinders
//==============================================================================

// The template parameters indicate the axes perpendicular to the axis of the
// cylinder.  offset1 and offset2 should correspond with i1 and i2,
// respectively.
template<int i1, int i2> double
axis_aligned_cylinder_evaluate(const double xyz[3], double offset1,
                               double offset2, double radius)
{
  const double xyz1 = xyz[i1] - offset1;
  const double xyz2 = xyz[i2] - offset2;
  return xyz1*xyz1 + xyz2*xyz2 - radius*radius;
}

// The first template parameter indicates which axis the cylinder is aligned to.
// The other two parameters indicate the other two axes.  offset1 and offset2
// should correspond with i2 and i3, respectively.
template<int i1, int i2, int i3> double
axis_aligned_cylinder_distance(const double xyz[3], const double uvw[3],
     bool coincident, double offset1, double offset2, double radius)
{
  const double a = 1.0 - uvw[i1]*uvw[i1];  // u^2 + v^2
  if (a == 0.0) return INFTY;

  const double xyz2 = xyz[i2] - offset1;
  const double xyz3 = xyz[i3] - offset2;
  const double k = xyz2 * uvw[i2] + xyz3 * uvw[i3];
  const double c = xyz2*xyz2 + xyz3*xyz3 - radius*radius;
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
template<int i1, int i2, int i3> void
axis_aligned_cylinder_normal(const double xyz[3], double uvw[3], double offset1,
                             double offset2)
{
  uvw[i2] = 2.0 * (xyz[i2] - offset1);
  uvw[i3] = 2.0 * (xyz[i3] - offset2);
  uvw[i1] = 0.0;
}

//==============================================================================
// SurfaceXCylinder implementation
//==============================================================================

SurfaceXCylinder::SurfaceXCylinder(pugi::xml_node surf_node)
  : Surface(surf_node)
{
  read_coeffs(surf_node, id, y0, z0, r);
}

inline double SurfaceXCylinder::evaluate(const double xyz[3]) const
{
  return axis_aligned_cylinder_evaluate<1, 2>(xyz, y0, z0, r);
}

inline double SurfaceXCylinder::distance(const double xyz[3],
     const double uvw[3], bool coincident) const
{
  return axis_aligned_cylinder_distance<0, 1, 2>(xyz, uvw, coincident, y0, z0,
                                                 r);
}

inline void SurfaceXCylinder::normal(const double xyz[3], double uvw[3]) const
{
  axis_aligned_cylinder_normal<0, 1, 2>(xyz, uvw, y0, z0);
}


void SurfaceXCylinder::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "x-cylinder");
  std::array<double, 3> coeffs {{y0, z0, r}};
  write_double_1D(group_id, "coefficients", coeffs);
}

//==============================================================================
// SurfaceYCylinder implementation
//==============================================================================

SurfaceYCylinder::SurfaceYCylinder(pugi::xml_node surf_node)
  : Surface(surf_node)
{
  read_coeffs(surf_node, id, x0, z0, r);
}

inline double SurfaceYCylinder::evaluate(const double xyz[3]) const
{
  return axis_aligned_cylinder_evaluate<0, 2>(xyz, x0, z0, r);
}

inline double SurfaceYCylinder::distance(const double xyz[3],
     const double uvw[3], bool coincident) const
{
  return axis_aligned_cylinder_distance<1, 0, 2>(xyz, uvw, coincident, x0, z0,
                                                 r);
}

inline void SurfaceYCylinder::normal(const double xyz[3], double uvw[3]) const
{
  axis_aligned_cylinder_normal<1, 0, 2>(xyz, uvw, x0, z0);
}

void SurfaceYCylinder::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "y-cylinder");
  std::array<double, 3> coeffs {{x0, z0, r}};
  write_double_1D(group_id, "coefficients", coeffs);
}

//==============================================================================
// SurfaceZCylinder implementation
//==============================================================================

SurfaceZCylinder::SurfaceZCylinder(pugi::xml_node surf_node)
  : Surface(surf_node)
{
  read_coeffs(surf_node, id, x0, y0, r);
}

inline double SurfaceZCylinder::evaluate(const double xyz[3]) const
{
  return axis_aligned_cylinder_evaluate<0, 1>(xyz, x0, y0, r);
}

inline double SurfaceZCylinder::distance(const double xyz[3],
     const double uvw[3], bool coincident) const
{
  return axis_aligned_cylinder_distance<2, 0, 1>(xyz, uvw, coincident, x0, y0,
                                                 r);
}

inline void SurfaceZCylinder::normal(const double xyz[3], double uvw[3]) const
{
  axis_aligned_cylinder_normal<2, 0, 1>(xyz, uvw, x0, y0);
}

void SurfaceZCylinder::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "z-cylinder");
  std::array<double, 3> coeffs {{x0, y0, r}};
  write_double_1D(group_id, "coefficients", coeffs);
}

//==============================================================================
// SurfaceSphere implementation
//==============================================================================

SurfaceSphere::SurfaceSphere(pugi::xml_node surf_node)
  : Surface(surf_node)
{
  read_coeffs(surf_node, id, x0, y0, z0, r);
}

double SurfaceSphere::evaluate(const double xyz[3]) const
{
  const double x = xyz[0] - x0;
  const double y = xyz[1] - y0;
  const double z = xyz[2] - z0;
  return x*x + y*y + z*z - r*r;
}

double SurfaceSphere::distance(const double xyz[3], const double uvw[3],
                               bool coincident) const
{
  const double x = xyz[0] - x0;
  const double y = xyz[1] - y0;
  const double z = xyz[2] - z0;
  const double k = x*uvw[0] + y*uvw[1] + z*uvw[2];
  const double c = x*x + y*y + z*z - r*r;
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

inline void SurfaceSphere::normal(const double xyz[3], double uvw[3]) const
{
  uvw[0] = 2.0 * (xyz[0] - x0);
  uvw[1] = 2.0 * (xyz[1] - y0);
  uvw[2] = 2.0 * (xyz[2] - z0);
}

void SurfaceSphere::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "sphere");
  std::array<double, 4> coeffs {{x0, y0, z0, r}};
  write_double_1D(group_id, "coefficients", coeffs);
}

//==============================================================================
// Generic functions for x-, y-, and z-, cones
//==============================================================================

// The first template parameter indicates which axis the cone is aligned to.
// The other two parameters indicate the other two axes.  offset1, offset2,
// and offset3 should correspond with i1, i2, and i3, respectively.
template<int i1, int i2, int i3> double
axis_aligned_cone_evaluate(const double xyz[3], double offset1,
                           double offset2, double offset3, double radius_sq)
{
  const double xyz1 = xyz[i1] - offset1;
  const double xyz2 = xyz[i2] - offset2;
  const double xyz3 = xyz[i3] - offset3;
  return xyz2*xyz2 + xyz3*xyz3 - radius_sq*xyz1*xyz1;
}

// The first template parameter indicates which axis the cone is aligned to.
// The other two parameters indicate the other two axes.  offset1, offset2,
// and offset3 should correspond with i1, i2, and i3, respectively.
template<int i1, int i2, int i3> double
axis_aligned_cone_distance(const double xyz[3], const double uvw[3],
     bool coincident, double offset1, double offset2, double offset3,
     double radius_sq)
{
  const double xyz1 = xyz[i1] - offset1;
  const double xyz2 = xyz[i2] - offset2;
  const double xyz3 = xyz[i3] - offset3;
  const double a = uvw[i2]*uvw[i2] + uvw[i3]*uvw[i3]
                   - radius_sq*uvw[i1]*uvw[i1];
  const double k = xyz2*uvw[i2] + xyz3*uvw[i3] - radius_sq*xyz1*uvw[i1];
  const double c = xyz2*xyz2 + xyz3*xyz3 - radius_sq*xyz1*xyz1;
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
template<int i1, int i2, int i3> void
axis_aligned_cone_normal(const double xyz[3], double uvw[3], double offset1,
                         double offset2, double offset3, double radius_sq)
{
  uvw[i1] = -2.0 * radius_sq * (xyz[i1] - offset1);
  uvw[i2] = 2.0 * (xyz[i2] - offset2);
  uvw[i3] = 2.0 * (xyz[i3] - offset3);
}

//==============================================================================
// SurfaceXCone implementation
//==============================================================================

SurfaceXCone::SurfaceXCone(pugi::xml_node surf_node)
  : Surface(surf_node)
{
  read_coeffs(surf_node, id, x0, y0, z0, r_sq);
}

inline double SurfaceXCone::evaluate(const double xyz[3]) const
{
  return axis_aligned_cone_evaluate<0, 1, 2>(xyz, x0, y0, z0, r_sq);
}

inline double SurfaceXCone::distance(const double xyz[3],
     const double uvw[3], bool coincident) const
{
  return axis_aligned_cone_distance<0, 1, 2>(xyz, uvw, coincident, x0, y0, z0,
                                             r_sq);
}

inline void SurfaceXCone::normal(const double xyz[3], double uvw[3]) const
{
  axis_aligned_cone_normal<0, 1, 2>(xyz, uvw, x0, y0, z0, r_sq);
}

void SurfaceXCone::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "x-cone");
  std::array<double, 4> coeffs {{x0, y0, z0, r_sq}};
  write_double_1D(group_id, "coefficients", coeffs);
}

//==============================================================================
// SurfaceYCone implementation
//==============================================================================

SurfaceYCone::SurfaceYCone(pugi::xml_node surf_node)
  : Surface(surf_node)
{
  read_coeffs(surf_node, id, x0, y0, z0, r_sq);
}

inline double SurfaceYCone::evaluate(const double xyz[3]) const
{
  return axis_aligned_cone_evaluate<1, 0, 2>(xyz, y0, x0, z0, r_sq);
}

inline double SurfaceYCone::distance(const double xyz[3],
     const double uvw[3], bool coincident) const
{
  return axis_aligned_cone_distance<1, 0, 2>(xyz, uvw, coincident, y0, x0, z0,
                                             r_sq);
}

inline void SurfaceYCone::normal(const double xyz[3], double uvw[3]) const
{
  axis_aligned_cone_normal<1, 0, 2>(xyz, uvw, y0, x0, z0, r_sq);
}

void SurfaceYCone::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "y-cone");
  std::array<double, 4> coeffs {{x0, y0, z0, r_sq}};
  write_double_1D(group_id, "coefficients", coeffs);
}

//==============================================================================
// SurfaceZCone implementation
//==============================================================================

SurfaceZCone::SurfaceZCone(pugi::xml_node surf_node)
  : Surface(surf_node)
{
  read_coeffs(surf_node, id, x0, y0, z0, r_sq);
}

inline double SurfaceZCone::evaluate(const double xyz[3]) const
{
  return axis_aligned_cone_evaluate<2, 0, 1>(xyz, z0, x0, y0, r_sq);
}

inline double SurfaceZCone::distance(const double xyz[3],
     const double uvw[3], bool coincident) const
{
  return axis_aligned_cone_distance<2, 0, 1>(xyz, uvw, coincident, z0, x0, y0,
                                             r_sq);
}

inline void SurfaceZCone::normal(const double xyz[3], double uvw[3]) const
{
  axis_aligned_cone_normal<2, 0, 1>(xyz, uvw, z0, x0, y0, r_sq);
}

void SurfaceZCone::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "z-cone");
  std::array<double, 4> coeffs {{x0, y0, z0, r_sq}};
  write_double_1D(group_id, "coefficients", coeffs);
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
SurfaceQuadric::evaluate(const double xyz[3]) const
{
  const double &x = xyz[0];
  const double &y = xyz[1];
  const double &z = xyz[2];
  return x*(A*x + D*y + G) +
         y*(B*y + E*z + H) +
         z*(C*z + F*x + J) + K;
}

double
SurfaceQuadric::distance(const double xyz[3],
     const double uvw[3], bool coincident) const
{
  const double &x = xyz[0];
  const double &y = xyz[1];
  const double &z = xyz[2];
  const double &u = uvw[0];
  const double &v = uvw[1];
  const double &w = uvw[2];

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

void
SurfaceQuadric::normal(const double xyz[3], double uvw[3]) const
{
  const double &x = xyz[0];
  const double &y = xyz[1];
  const double &z = xyz[2];
  uvw[0] = 2.0*A*x + D*y + F*z + G;
  uvw[1] = 2.0*B*y + D*x + E*z + H;
  uvw[2] = 2.0*C*z + E*y + F*x + J;
}

void SurfaceQuadric::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "quadric");
  std::array<double, 10> coeffs {{A, B, C, D, E, F, G, H, J, K}};
  write_double_1D(group_id, "coefficients", coeffs);
}

//==============================================================================

extern "C" void
read_surfaces(pugi::xml_node *node)
{
  // Count the number of surfaces.
  for (pugi::xml_node surf_node: node->children("surface")) {n_surfaces++;}
  if (n_surfaces == 0) {
    fatal_error("No surfaces found in geometry.xml!");
  }

  // Allocate the array of Surface pointers.
  surfaces_c = new Surface* [n_surfaces];

  // Loop over XML surface elements and populate the array.
  {
    pugi::xml_node surf_node;
    int i_surf;
    for (surf_node = node->child("surface"), i_surf = 0; surf_node;
         surf_node = surf_node.next_sibling("surface"), i_surf++) {
      std::string surf_type = get_node_value(surf_node, "type");

      if (surf_type == "x-plane") {
        surfaces_c[i_surf] = new SurfaceXPlane(surf_node);

      } else if (surf_type == "y-plane") {
        surfaces_c[i_surf] = new SurfaceYPlane(surf_node);

      } else if (surf_type == "z-plane") {
        surfaces_c[i_surf] = new SurfaceZPlane(surf_node);

      } else if (surf_type == "plane") {
        surfaces_c[i_surf] = new SurfacePlane(surf_node);

      } else if (surf_type == "x-cylinder") {
        surfaces_c[i_surf] = new SurfaceXCylinder(surf_node);

      } else if (surf_type == "y-cylinder") {
        surfaces_c[i_surf] = new SurfaceYCylinder(surf_node);

      } else if (surf_type == "z-cylinder") {
        surfaces_c[i_surf] = new SurfaceZCylinder(surf_node);

      } else if (surf_type == "sphere") {
        surfaces_c[i_surf] = new SurfaceSphere(surf_node);

      } else if (surf_type == "x-cone") {
        surfaces_c[i_surf] = new SurfaceXCone(surf_node);

      } else if (surf_type == "y-cone") {
        surfaces_c[i_surf] = new SurfaceYCone(surf_node);

      } else if (surf_type == "z-cone") {
        surfaces_c[i_surf] = new SurfaceZCone(surf_node);

      } else if (surf_type == "quadric") {
        surfaces_c[i_surf] = new SurfaceQuadric(surf_node);

      } else {
        std::stringstream err_msg;
        err_msg << "Invalid surface type, \"" << surf_type << "\"";
        fatal_error(err_msg);
      }
    }
  }

  // Fill the surface dictionary.
  for (int i_surf = 0; i_surf < n_surfaces; i_surf++) {
    int id = surfaces_c[i_surf]->id;
    auto in_dict = surface_dict.find(id);
    if (in_dict == surface_dict.end()) {
      surface_dict[id] = i_surf;
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
    if (surfaces_c[i_surf]->bc == BC_PERIODIC) {
      // Downcast to the PeriodicSurface type.
      Surface *surf_base = surfaces_c[i_surf];
      PeriodicSurface *surf = dynamic_cast<PeriodicSurface *>(surf_base);

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
    if (surfaces_c[i_surf]->bc == BC_PERIODIC) {
      // Downcast to the PeriodicSurface type.
      Surface *surf_base = surfaces_c[i_surf];
      PeriodicSurface *surf = dynamic_cast<PeriodicSurface *>(surf_base);

      // Also try downcasting to the SurfacePlane type (which must be handled
      // differently).
      SurfacePlane *surf_p = dynamic_cast<SurfacePlane *>(surf);

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
          surf->i_periodic = surface_dict[surf->i_periodic];
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
          surf->i_periodic = surface_dict[surf->i_periodic];
        }
      }

      // Make sure the opposite surface is also periodic.
      if (surfaces_c[surf->i_periodic]->bc != BC_PERIODIC) {
        std::stringstream err_msg;
        err_msg << "Could not find matching surface for periodic boundary "
                   "condition on surface " << surf->id;
        fatal_error(err_msg);
      }
    }
  }
}

} // namespace openmc
