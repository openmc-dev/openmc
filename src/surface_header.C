#include <cstring>  // For strcmp
#include <iostream>
#include <limits>  // For numeric_limits
#include <math.h>  // For fabs
#include "pugixml/pugixml.hpp"

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
// SURFACE type defines a first- or second-order surface that can be used to
// construct closed volumes (cells)
//==============================================================================

class Surface {
  int id;                    // Unique ID
  int neighbor_pos[],        // List of cells on positive side
      neighbor_neg[];        // List of cells on negative side
  int bc;                    // Boundary condition
  //TODO: swith that zero to a NONE constant.
  int i_periodic = 0;        // Index of corresponding periodic surface
  char name[104];            // User-defined name

public:
  bool sense(const double xyz[3], const double uvw[3]) const;
  void reflect(const double xyz[3], double uvw[3]) const;
  virtual double evaluate(const double xyz[3]) const = 0;
  virtual double distance(const double xyz[3], const double uvw[3],
                          bool coincident) const = 0;
  virtual void normal(const double xyz[3], double uvw[3]) const = 0;
};

bool
Surface::sense(const double xyz[3], const double uvw[3]) const {
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
Surface::reflect(const double xyz[3], double uvw[3]) const {
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
// Generic functions for X-, Y-, and Z-, planes
//==============================================================================

template<int i> double
axis_aligned_plane_evaluate(const double xyz[3], double offset) {
  return xyz[i] - offset;}

template<int i> double
axis_aligned_plane_distance(const double xyz[3], const double uvw[3],
                            bool coincident, double offset) {
  const double f = offset - xyz[i];
  if (coincident or fabs(f) < FP_COINCIDENT or uvw[i] == 0.0) return INFTY;
  const double d = f / uvw[i];
  if (d < 0.0) return INFTY;
  return d;
}

template<int i1, int i2, int i3> void
axis_aligned_plane_normal(const double xyz[3], double uvw[3]) {
  uvw[i1] = 1.0;
  uvw[i2] = 0.0;
  uvw[i3] = 0.0;
}

//==============================================================================
// SurfaceXPlane
//==============================================================================

class SurfaceXPlane : public Surface {
  double x0;
public:
  SurfaceXPlane(pugi::xml_node surf_node);
  double evaluate(const double xyz[3]) const;
  double distance(const double xyz[3], const double uvw[3],
                  const bool coincident) const;
  void normal(const double xyz[3], double uvw[3]) const;
};

SurfaceXPlane::SurfaceXPlane(pugi::xml_node surf_node) {
  const char *coeffs = surf_node.attribute("coeffs").value();
  int stat = sscanf(coeffs, "%lf", &x0);
  if (stat != 1) {
    std::cout << "Something went wrong reading surface coeffs!" << std::endl;
  }
}

inline double SurfaceXPlane::evaluate(const double xyz[3]) const {
  return axis_aligned_plane_evaluate<0>(xyz, x0);
}

inline double SurfaceXPlane::distance(const double xyz[3], const double uvw[3],
                                      bool coincident) const {
  return axis_aligned_plane_distance<0>(xyz, uvw, coincident, x0);
}

inline void SurfaceXPlane::normal(const double xyz[3], double uvw[3]) const {
  axis_aligned_plane_normal<0, 1, 2>(xyz, uvw);
}

//==============================================================================
// SurfaceYPlane
//==============================================================================

class SurfaceYPlane : public Surface {
  double y0;
public:
  SurfaceYPlane(pugi::xml_node surf_node);
  double evaluate(const double xyz[3]) const;
  double distance(const double xyz[3], const double uvw[3],
                  bool coincident) const;
  void normal(const double xyz[3], double uvw[3]) const;
};

SurfaceYPlane::SurfaceYPlane(pugi::xml_node surf_node) {
  const char *coeffs = surf_node.attribute("coeffs").value();
  int stat = sscanf(coeffs, "%lf", &y0);
  if (stat != 1) {
    std::cout << "Something went wrong reading surface coeffs!" << std::endl;
  }
}

inline double SurfaceYPlane::evaluate(const double xyz[3]) const {
  return axis_aligned_plane_evaluate<1>(xyz, y0);
}

inline double SurfaceYPlane::distance(const double xyz[3], const double uvw[3],
                                      bool coincident) const {
  return axis_aligned_plane_distance<1>(xyz, uvw, coincident, y0);
}

inline void SurfaceYPlane::normal(const double xyz[3], double uvw[3]) const {
  axis_aligned_plane_normal<1, 0, 2>(xyz, uvw);
}

//==============================================================================
// SurfaceZPlane
//==============================================================================

class SurfaceZPlane : public Surface {
  double z0;
public:
  SurfaceZPlane(pugi::xml_node surf_node);
  double evaluate(const double xyz[3]) const;
  double distance(const double xyz[3], const double uvw[3],
                  bool coincident) const;
  void normal(const double xyz[3], double uvw[3]) const;
};

SurfaceZPlane::SurfaceZPlane(pugi::xml_node surf_node) {
  const char *coeffs = surf_node.attribute("coeffs").value();
  int stat = sscanf(coeffs, "%lf", &z0);
  if (stat != 1) {
    std::cout << "Something went wrong reading surface coeffs!" << std::endl;
  }
}

inline double SurfaceZPlane::evaluate(const double xyz[3]) const {
  return axis_aligned_plane_evaluate<2>(xyz, z0);
}

inline double SurfaceZPlane::distance(const double xyz[3], const double uvw[3],
                                      bool coincident) const {
  return axis_aligned_plane_distance<2>(xyz, uvw, coincident, z0);
}

inline void SurfaceZPlane::normal(const double xyz[3], double uvw[3]) const {
  axis_aligned_plane_normal<2, 0, 1>(xyz, uvw);
}

//==============================================================================
// SurfacePlane
//==============================================================================

// TODO

//==============================================================================
// Generic functions for X-, Y-, and Z-, cylinders
//==============================================================================

template<int i1, int i2> double
axis_aligned_cylinder_evaluate(const double xyz[3], double offset1,
                               double offset2, double radius) {
  const double xyz1 = xyz[i1] - offset1;
  const double xyz2 = xyz[i2] - offset2;
  return xyz1*xyz1 + xyz2*xyz2 - radius*radius;
}

template<int i1, int i2, int i3> double
axis_aligned_cylinder_distance(const double xyz[3], const double uvw[3],
     double offset1, double offset2, double radius, bool coincident) {
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

template<int i1, int i2, int i3> void
axis_aligned_cylinder_normal(const double xyz[3], double uvw[3], double offset1,
                             double offset2) {
  uvw[i2] = 2.0 * (xyz[i2] - offset1);
  uvw[i3] = 2.0 * (xyz[i3] - offset2);
  uvw[i1] = 0.0;
}

//==============================================================================
// SurfaceXCylinder
//==============================================================================

class SurfaceXCylinder : public Surface {
  double y0, z0, r;
public:
  SurfaceXCylinder(pugi::xml_node surf_node);
  double evaluate(const double xyz[3]) const;
  double distance(const double xyz[3], const double uvw[3],
                  bool coincident) const;
  void normal(const double xyz[3], double uvw[3]) const;
};

SurfaceXCylinder::SurfaceXCylinder(pugi::xml_node surf_node) {
  const char *coeffs = surf_node.attribute("coeffs").value();
  int stat = sscanf(coeffs, "%lf %lf %lf", &y0, &z0, &r);
  if (stat != 3) {
    std::cout << "Something went wrong reading surface coeffs!" << std::endl;
  }
}

inline double
SurfaceXCylinder::evaluate(const double xyz[3]) const {
  return axis_aligned_cylinder_evaluate<1, 2>(xyz, y0, z0, r);
}

inline double SurfaceXCylinder::distance(const double xyz[3],
     const double uvw[3], bool coincident) const {
  return axis_aligned_cylinder_distance<0, 1, 2>(xyz, uvw, coincident, y0, z0,
                                                 r);
}

inline void SurfaceXCylinder::normal(const double xyz[3], double uvw[3]) const {
  axis_aligned_cylinder_normal<0, 1, 2>(xyz, uvw, y0, z0);
}

//==============================================================================
// SurfaceYCylinder
//==============================================================================

class SurfaceYCylinder : public Surface {
  double x0, z0, r;
public:
  SurfaceYCylinder(pugi::xml_node surf_node);
  double evaluate(const double xyz[3]) const;
  double distance(const double xyz[3], const double uvw[3],
                  bool coincident) const;
  void normal(const double xyz[3], double uvw[3]) const;
};

SurfaceYCylinder::SurfaceYCylinder(pugi::xml_node surf_node) {
  const char *coeffs = surf_node.attribute("coeffs").value();
  int stat = sscanf(coeffs, "%lf %lf %lf", &x0, &z0, &r);
  if (stat != 3) {
    std::cout << "Something went wrong reading surface coeffs!" << std::endl;
  }
}

inline double
SurfaceYCylinder::evaluate(const double xyz[3]) const {
  return axis_aligned_cylinder_evaluate<0, 2>(xyz, x0, z0, r);
}

inline double SurfaceYCylinder::distance(const double xyz[3],
     const double uvw[3], bool coincident) const {
  return axis_aligned_cylinder_distance<1, 0, 2>(xyz, uvw, coincident, x0, z0,
                                                 r);
}

inline void SurfaceYCylinder::normal(const double xyz[3], double uvw[3]) const {
  axis_aligned_cylinder_normal<1, 0, 2>(xyz, uvw, x0, z0);
}

//==============================================================================
// SurfaceZCylinder
//==============================================================================

class SurfaceZCylinder : public Surface {
  double x0, y0, r;
public:
  SurfaceZCylinder(pugi::xml_node surf_node);
  double evaluate(const double xyz[3]) const;
  double distance(const double xyz[3], const double uvw[3],
                  bool coincident) const;
  void normal(const double xyz[3], double uvw[3]) const;
};

SurfaceZCylinder::SurfaceZCylinder(pugi::xml_node surf_node) {
  const char *coeffs = surf_node.attribute("coeffs").value();
  int stat = sscanf(coeffs, "%lf %lf %lf", &x0, &y0, &r);
  if (stat != 3) {
    std::cout << "Something went wrong reading surface coeffs!" << std::endl;
  }
}

inline double
SurfaceZCylinder::evaluate(const double xyz[3]) const {
  return axis_aligned_cylinder_evaluate<0, 1>(xyz, x0, y0, r);
}

inline double SurfaceZCylinder::distance(const double xyz[3],
     const double uvw[3], bool coincident) const {
  return axis_aligned_cylinder_distance<2, 0, 1>(xyz, uvw, coincident, x0, y0,
                                                 r);
}

inline void SurfaceZCylinder::normal(const double xyz[3], double uvw[3]) const {
  axis_aligned_cylinder_normal<2, 0, 1>(xyz, uvw, x0, y0);
}

//==============================================================================
// SurfaceSphere
//==============================================================================

class SurfaceSphere : public Surface {
  double x0, y0, z0, r;
public:
  SurfaceSphere(pugi::xml_node surf_node);
  double evaluate(const double xyz[3]) const;
  double distance(const double xyz[3], const double uvw[3],
                  bool coincident) const;
  void normal(const double xyz[3], double uvw[3]) const;
};

SurfaceSphere::SurfaceSphere(pugi::xml_node surf_node) {
  const char *coeffs = surf_node.attribute("coeffs").value();
  int stat = sscanf(coeffs, "%lf %lf %lf %lf", &x0, &y0, &z0, &r);
  if (stat != 4) {
    std::cout << "Something went wrong reading surface coeffs!" << std::endl;
  }
}

double SurfaceSphere::evaluate(const double xyz[3]) const {
  const double x = xyz[0] - x0;
  const double y = xyz[1] - y0;
  const double z = xyz[2] - z0;
  return x*x + y*y + z*z - r*r;
}

double SurfaceSphere::distance(const double xyz[3], const double uvw[3],
                               bool coincident) const {
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

inline void SurfaceSphere::normal(const double xyz[3], double uvw[3]) const {
  uvw[0] = 2.0 * (xyz[0] - x0);
  uvw[1] = 2.0 * (xyz[1] - y0);
  uvw[2] = 2.0 * (xyz[2] - z0);
}

//==============================================================================
// SurfaceXCone
//==============================================================================

// TODO

//==============================================================================
// SurfaceYCone
//==============================================================================

// TODO

//==============================================================================
// SurfaceZCone
//==============================================================================

// TODO

//==============================================================================
// SurfaceQuadric
//==============================================================================

// TODO

//==============================================================================

extern "C" void
read_surfaces(pugi::xml_node *node) {
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
      if (surf_node.attribute("type")) {
        const pugi::char_t *surf_type = surf_node.attribute("type").value();

        if (strcmp(surf_type, "x-plane") == 0) {
          surfaces_c[i_surf] = new SurfaceXPlane(surf_node);

        } else if (strcmp(surf_type, "y-plane") == 0) {
          surfaces_c[i_surf] = new SurfaceYPlane(surf_node);

        } else if (strcmp(surf_type, "z-plane") == 0) {
          surfaces_c[i_surf] = new SurfaceZPlane(surf_node);

        } else if (strcmp(surf_type, "x-cylinder") == 0) {
          surfaces_c[i_surf] = new SurfaceXCylinder(surf_node);

        } else if (strcmp(surf_type, "y-cylinder") == 0) {
          surfaces_c[i_surf] = new SurfaceYCylinder(surf_node);

        } else if (strcmp(surf_type, "z-cylinder") == 0) {
          surfaces_c[i_surf] = new SurfaceZCylinder(surf_node);

        } else if (strcmp(surf_type, "sphere") == 0) {
          surfaces_c[i_surf] = new SurfaceSphere(surf_node);

        } else {
          std::cout << "Call error or handle uppercase here!" << std::endl;
          std::cout << surf_type << std::endl;
        }
      }
    }
  }
}

//==============================================================================

extern "C" bool
surface_sense(int surf_ind, double xyz[3], double uvw[3]) {
  return surfaces_c[surf_ind]->sense(xyz, uvw);
}

extern "C" double
surface_distance(int surf_ind, double xyz[3], double uvw[3], bool coincident) {
  return surfaces_c[surf_ind]->distance(xyz, uvw, coincident);
}
