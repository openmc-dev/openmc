#include <cstring>  // For strcmp
#include <limits>  // For numeric_limits
#include <math.h>  // For fabs

#include "pugixml/pugixml.hpp"

// DEBUGGING
#include <iostream>

//==============================================================================
// Constants
//==============================================================================

const double FP_COINCIDENT = 1e-12;
const double INFTY = std::numeric_limits<double>::max();

//==============================================================================
// Global array of surfaces
//==============================================================================

class Surface;
Surface **surfaces_c;

//==============================================================================
// Helper functions for reading the "coeffs" node of an XML surface element
//==============================================================================

const char* get_coeff_str(pugi::xml_node surf_node)
{
  if (surf_node.attribute("coeffs")) {
    return surf_node.attribute("coeffs").value();
  } else if (surf_node.child("coeffs")) {
    return surf_node.child_value("coeffs");
  } else {
    std::cout << "ERROR: Found a surface with no coefficients" << std::endl;
    return NULL;
  }
}

void read_coeffs(pugi::xml_node surf_node, double &c1)
{
  const char *coeffs = get_coeff_str(surf_node);
  int stat = sscanf(coeffs, "%lf", &c1);
  if (stat != 1) {
    std::cout << "Something went wrong reading surface coeffs!" << std::endl;
  }
}

void read_coeffs(pugi::xml_node surf_node, double &c1, double &c2, double &c3)
{
  const char *coeffs = get_coeff_str(surf_node);
  int stat = sscanf(coeffs, "%lf %lf %lf", &c1, &c2, &c3);
  if (stat != 3) {
    std::cout << "Something went wrong reading surface coeffs!" << std::endl;
  }
}

void read_coeffs(pugi::xml_node surf_node, double &c1, double &c2, double &c3,
                 double &c4)
{
  const char *coeffs = get_coeff_str(surf_node);
  int stat = sscanf(coeffs, "%lf %lf %lf %lf", &c1, &c2, &c3, &c4);
  if (stat != 4) {
    std::cout << "Something went wrong reading surface coeffs!" << std::endl;
  }
}

void read_coeffs(pugi::xml_node surf_node, double &c1, double &c2, double &c3,
                 double &c4, double &c5, double &c6, double &c7, double &c8,
                 double &c9, double &c10)
{
  const char *coeffs = get_coeff_str(surf_node);
  int stat = sscanf(coeffs, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                    &c1, &c2, &c3, &c4, &c5, &c6, &c7, &c8, &c9, &c10);
  if (stat != 10) {
    std::cout << "Something went wrong reading surface coeffs!" << std::endl;
  }
}

//==============================================================================
//! A geometry primitive used to define regions of 3D space.
//==============================================================================

class Surface
{
public:
  int id;                    //!< Unique ID
  int neighbor_pos[],        //!< List of cells on positive side
      neighbor_neg[];        //!< List of cells on negative side
  int bc;                    //!< Boundary condition
  //TODO: switch that zero to a NONE constant.
  int i_periodic = 0;        //!< Index of corresponding periodic surface
  char name[104];            //!< User-defined name

  //! Determine which side of a surface a point lies on.
  //! @param xyz[3] The 3D Cartesian coordinate of a point.
  //! @param uvw[3] A direction used to "break ties" and pick a sense when the
  //!   point is very close to the surface.
  //! @return True if the point is on the "positive" side of the surface and
  //!   false otherwise.
  bool sense(const double xyz[3], const double uvw[3]) const;

  //! Determine the direction of a ray reflected from the surface.
  //! @param xyz[3] The point at which the ray is incident.
  //! @param uvw[3] A direction.  This is both an input and an output parameter.
  //!   It specifies the icident direction on input and the reflected direction
  //!   on output.
  void reflect(const double xyz[3], double uvw[3]) const;

  //! Evaluate the equation describing the surface.
  //!
  //! Surfaces can be described by some function f(x, y, z) = 0.  This member
  //! function evaluates that mathematical function.
  //! @param xyz[3] A 3D Cartesian coordinate.
  virtual double evaluate(const double xyz[3]) const = 0;

  //! Compute the distance between a point and the surface along a ray.
  //! @param xyz[3] A 3D Cartesian coordinate.
  //! @param uvw[3] The direction of the ray.
  //! @param coincident A hint to the code that the given point should lie
  //!   exactly on the surface.
  virtual double distance(const double xyz[3], const double uvw[3],
                          bool coincident) const = 0;

  //! Compute the local outward normal direction of the surface.
  //! @param xyz[3] A 3D Cartesian coordinate.
  //! @param uvw[3] This output argument provides the normal.
  virtual void normal(const double xyz[3], double uvw[3]) const = 0;
};

bool
Surface::sense(const double xyz[3], const double uvw[3]) const
{
  // Evaluate the surface equation at the particle's coordinates to determine
  // which side the particle is on.
  const double f = evaluate(xyz);

  // Check which side of surface the point is on.
  if (fabs(f) < FP_COINCIDENT) {
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
  if (coincident or fabs(f) < FP_COINCIDENT or uvw[i] == 0.0) return INFTY;
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
// SurfaceXPlane
//! A plane perpendicular to the x-axis.
//
//! The plane is described by the equation \f$x - x_0 = 0\f$
//==============================================================================

class SurfaceXPlane : public Surface
{
  double x0;
public:
  SurfaceXPlane(pugi::xml_node surf_node);
  double evaluate(const double xyz[3]) const;
  double distance(const double xyz[3], const double uvw[3],
                  const bool coincident) const;
  void normal(const double xyz[3], double uvw[3]) const;
};

SurfaceXPlane::SurfaceXPlane(pugi::xml_node surf_node)
{
  read_coeffs(surf_node, x0);
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

//==============================================================================
// SurfaceYPlane
//! A plane perpendicular to the y-axis.
//
//! The plane is described by the equation \f$y - y_0 = 0\f$
//==============================================================================

class SurfaceYPlane : public Surface
{
  double y0;
public:
  SurfaceYPlane(pugi::xml_node surf_node);
  double evaluate(const double xyz[3]) const;
  double distance(const double xyz[3], const double uvw[3],
                  bool coincident) const;
  void normal(const double xyz[3], double uvw[3]) const;
};

SurfaceYPlane::SurfaceYPlane(pugi::xml_node surf_node)
{
  read_coeffs(surf_node, y0);
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

//==============================================================================
// SurfaceZPlane
//! A plane perpendicular to the z-axis.
//
//! The plane is described by the equation \f$z - z_0 = 0\f$
//==============================================================================

class SurfaceZPlane : public Surface
{
  double z0;
public:
  SurfaceZPlane(pugi::xml_node surf_node);
  double evaluate(const double xyz[3]) const;
  double distance(const double xyz[3], const double uvw[3],
                  bool coincident) const;
  void normal(const double xyz[3], double uvw[3]) const;
};

SurfaceZPlane::SurfaceZPlane(pugi::xml_node surf_node)
{
  read_coeffs(surf_node, z0);
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

//==============================================================================
// SurfacePlane
//! A general plane.
//
//! The plane is described by the equation \f$A x + B y + C z - D = 0\f$
//==============================================================================

class SurfacePlane : public Surface
{
  double A, B, C, D;
public:
  SurfacePlane(pugi::xml_node surf_node);
  double evaluate(const double xyz[3]) const;
  double distance(const double xyz[3], const double uvw[3],
                  const bool coincident) const;
  void normal(const double xyz[3], double uvw[3]) const;
};

SurfacePlane::SurfacePlane(pugi::xml_node surf_node)
{
  read_coeffs(surf_node, A, B, C, D);
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
  if (coincident or fabs(f) < FP_COINCIDENT or projection == 0.0) {
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

  } else if (coincident or fabs(c) < FP_COINCIDENT) {
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
// SurfaceXCylinder
//! A cylinder aligned along the x-axis.
//
//! The cylinder is described by the equation
//! \f$(y - y_0)^2 + (z - z_0)^2 - R^2 = 0\f$
//==============================================================================

class SurfaceXCylinder : public Surface
{
  double y0, z0, r;
public:
  SurfaceXCylinder(pugi::xml_node surf_node);
  double evaluate(const double xyz[3]) const;
  double distance(const double xyz[3], const double uvw[3],
                  bool coincident) const;
  void normal(const double xyz[3], double uvw[3]) const;
};

SurfaceXCylinder::SurfaceXCylinder(pugi::xml_node surf_node)
{
  read_coeffs(surf_node, y0, z0, r);
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

//==============================================================================
// SurfaceYCylinder
//! A cylinder aligned along the y-axis.
//
//! The cylinder is described by the equation
//! \f$(x - x_0)^2 + (z - z_0)^2 - R^2 = 0\f$
//==============================================================================

class SurfaceYCylinder : public Surface
{
  double x0, z0, r;
public:
  SurfaceYCylinder(pugi::xml_node surf_node);
  double evaluate(const double xyz[3]) const;
  double distance(const double xyz[3], const double uvw[3],
                  bool coincident) const;
  void normal(const double xyz[3], double uvw[3]) const;
};

SurfaceYCylinder::SurfaceYCylinder(pugi::xml_node surf_node)
{
  read_coeffs(surf_node, x0, z0, r);
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

//==============================================================================
// SurfaceZCylinder
//! A cylinder aligned along the z-axis.
//
//! The cylinder is described by the equation
//! \f$(x - x_0)^2 + (y - y_0)^2 - R^2 = 0\f$
//==============================================================================

class SurfaceZCylinder : public Surface
{
  double x0, y0, r;
public:
  SurfaceZCylinder(pugi::xml_node surf_node);
  double evaluate(const double xyz[3]) const;
  double distance(const double xyz[3], const double uvw[3],
                  bool coincident) const;
  void normal(const double xyz[3], double uvw[3]) const;
};

SurfaceZCylinder::SurfaceZCylinder(pugi::xml_node surf_node)
{
  read_coeffs(surf_node, x0, y0, r);
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

//==============================================================================
// SurfaceSphere
//! A sphere.
//
//! The cylinder is described by the equation
//! \f$(x - x_0)^2 + (y - y_0)^2 + (z - z_0)^2 - R^2 = 0\f$
//==============================================================================

class SurfaceSphere : public Surface
{
  double x0, y0, z0, r;
public:
  SurfaceSphere(pugi::xml_node surf_node);
  double evaluate(const double xyz[3]) const;
  double distance(const double xyz[3], const double uvw[3],
                  bool coincident) const;
  void normal(const double xyz[3], double uvw[3]) const;
};

SurfaceSphere::SurfaceSphere(pugi::xml_node surf_node)
{
  read_coeffs(surf_node, x0, y0, z0, r);
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

  } else if (coincident or fabs(c) < FP_COINCIDENT) {
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

  } else if (coincident or fabs(c) < FP_COINCIDENT) {
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
// SurfaceXCone
//! A cone aligned along the x-axis.
//
//! The cylinder is described by the equation
//! \f$(y - y_0)^2 + (z - z_0)^2 - R^2 (x - x_0)^2 = 0\f$
//==============================================================================

class SurfaceXCone : public Surface
{
  double x0, y0, z0, r_sq;
public:
  SurfaceXCone(pugi::xml_node surf_node);
  double evaluate(const double xyz[3]) const;
  double distance(const double xyz[3], const double uvw[3],
                  bool coincident) const;
  void normal(const double xyz[3], double uvw[3]) const;
};

SurfaceXCone::SurfaceXCone(pugi::xml_node surf_node)
{
  read_coeffs(surf_node, x0, y0, z0, r_sq);
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

//==============================================================================
// SurfaceYCone
//! A cone aligned along the y-axis.
//
//! The cylinder is described by the equation
//! \f$(x - x_0)^2 + (z - z_0)^2 - R^2 (y - y_0)^2 = 0\f$
//==============================================================================

class SurfaceYCone : public Surface
{
  double x0, y0, z0, r_sq;
public:
  SurfaceYCone(pugi::xml_node surf_node);
  double evaluate(const double xyz[3]) const;
  double distance(const double xyz[3], const double uvw[3],
                  bool coincident) const;
  void normal(const double xyz[3], double uvw[3]) const;
};

SurfaceYCone::SurfaceYCone(pugi::xml_node surf_node)
{
  read_coeffs(surf_node, x0, y0, z0, r_sq);
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

//==============================================================================
// SurfaceZCone
//! A cone aligned along the z-axis.
//
//! The cylinder is described by the equation
//! \f$(x - x_0)^2 + (y - y_0)^2 - R^2 (z - z_0)^2 = 0\f$
//==============================================================================

class SurfaceZCone : public Surface
{
  double x0, y0, z0, r_sq;
public:
  SurfaceZCone(pugi::xml_node surf_node);
  double evaluate(const double xyz[3]) const;
  double distance(const double xyz[3], const double uvw[3],
                  bool coincident) const;
  void normal(const double xyz[3], double uvw[3]) const;
};

SurfaceZCone::SurfaceZCone(pugi::xml_node surf_node)
{
  read_coeffs(surf_node, x0, y0, z0, r_sq);
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

//==============================================================================
// SurfaceQuadric
//! A general surface described by a quadratic equation.
//
//! \f$A x^2 + B y^2 + C z^2 + D x y + E y z + F x z + G x + H y + J z + K = 0\f$
//==============================================================================

class SurfaceQuadric : public Surface
{
  // Ax^2 + By^2 + Cz^2 + Dxy + Eyz + Fxz + Gx + Hy + Jz + K = 0
  double A, B, C, D, E, F, G, H, J, K;
public:
  SurfaceQuadric(pugi::xml_node surf_node);
  double evaluate(const double xyz[3]) const;
  double distance(const double xyz[3], const double uvw[3],
                  bool coincident) const;
  void normal(const double xyz[3], double uvw[3]) const;
};

SurfaceQuadric::SurfaceQuadric(pugi::xml_node surf_node)
{
  read_coeffs(surf_node, A, B, C, D, E, F, G, H, J, K);
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

  } else if (coincident or fabs(c) < FP_COINCIDENT) {
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

//==============================================================================

extern "C" void
read_surfaces(pugi::xml_node *node)
{
  // Count the number of surfaces.
  int n_surfaces = 0;
  for (pugi::xml_node surf_node = node->child("surface"); surf_node;
       surf_node = surf_node.next_sibling("surface")) {
    n_surfaces += 1;
  }

  // Allocate the array of Surface pointers.
  surfaces_c = new Surface* [n_surfaces];

  // Loop over XML surface elements and populate the array.
  {
    pugi::xml_node surf_node;
    int i_surf;
    for (surf_node = node->child("surface"), i_surf = 0; surf_node;
         surf_node = surf_node.next_sibling("surface"), i_surf++) {
      const pugi::char_t *surf_type;
      if (surf_node.attribute("type")) {
        surf_type = surf_node.attribute("type").value();
      } else if (surf_node.child("type")) {
        surf_type = surf_node.child_value("type");
      } else {
        std::cout << "ERROR: Found a surface with no type attribute/node"
             << std::endl;
      }

      if (strcmp(surf_type, "x-plane") == 0) {
        surfaces_c[i_surf] = new SurfaceXPlane(surf_node);

      } else if (strcmp(surf_type, "y-plane") == 0) {
        surfaces_c[i_surf] = new SurfaceYPlane(surf_node);

      } else if (strcmp(surf_type, "z-plane") == 0) {
        surfaces_c[i_surf] = new SurfaceZPlane(surf_node);

      } else if (strcmp(surf_type, "plane") == 0) {
        surfaces_c[i_surf] = new SurfacePlane(surf_node);

      } else if (strcmp(surf_type, "x-cylinder") == 0) {
        surfaces_c[i_surf] = new SurfaceXCylinder(surf_node);

      } else if (strcmp(surf_type, "y-cylinder") == 0) {
        surfaces_c[i_surf] = new SurfaceYCylinder(surf_node);

      } else if (strcmp(surf_type, "z-cylinder") == 0) {
        surfaces_c[i_surf] = new SurfaceZCylinder(surf_node);

      } else if (strcmp(surf_type, "sphere") == 0) {
        surfaces_c[i_surf] = new SurfaceSphere(surf_node);

      } else if (strcmp(surf_type, "x-cone") == 0) {
        surfaces_c[i_surf] = new SurfaceXCone(surf_node);

      } else if (strcmp(surf_type, "y-cone") == 0) {
        surfaces_c[i_surf] = new SurfaceYCone(surf_node);

      } else if (strcmp(surf_type, "z-cone") == 0) {
        surfaces_c[i_surf] = new SurfaceZCone(surf_node);

      } else if (strcmp(surf_type, "quadric") == 0) {
        surfaces_c[i_surf] = new SurfaceQuadric(surf_node);

      } else {
        std::cout << "Call error or handle uppercase here!" << std::endl;
        std::cout << surf_type << std::endl;
      }
    }
  }
}

//==============================================================================

extern "C" bool
surface_sense(int surf_ind, double xyz[3], double uvw[3])
{
  return surfaces_c[surf_ind]->sense(xyz, uvw);
}

extern "C" double
surface_distance(int surf_ind, double xyz[3], double uvw[3], bool coincident)
{
  return surfaces_c[surf_ind]->distance(xyz, uvw, coincident);
}

extern "C" void
surface_normal(int surf_ind, double xyz[3], double uvw[3])
{
  return surfaces_c[surf_ind]->normal(xyz, uvw);
}
