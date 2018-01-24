#ifndef SURFACE_H
#define SURFACE_H

#include <map>
#include <limits>  // For numeric_limits

#include "hdf5.h"
#include "pugixml/pugixml.hpp"

//==============================================================================
// Module constants
//==============================================================================

extern "C" const int BC_TRANSMIT{0};
extern "C" const int BC_VACUUM{1};
extern "C" const int BC_REFLECT{2};
extern "C" const int BC_PERIODIC{3};

//==============================================================================
// Constants that should eventually be moved out of this file
//==============================================================================

extern "C" double FP_COINCIDENT;
const double INFTY{std::numeric_limits<double>::max()};
const int C_NONE{-1};

//==============================================================================
// Global variables
//==============================================================================

class Surface;
Surface **surfaces_c;

std::map<int, int> surface_dict;

//==============================================================================
//! Coordinates for an axis-aligned cube that bounds a geometric object.
//==============================================================================

struct BoundingBox
{
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  double zmin;
  double zmax;
};

//==============================================================================
//! A geometry primitive used to define regions of 3D space.
//==============================================================================

class Surface
{
public:
  int id;                    //!< Unique ID
  //int neighbor_pos[],        //!< List of cells on positive side
  //    neighbor_neg[];        //!< List of cells on negative side
  int bc;                    //!< Boundary condition
  std::string name{""};      //!< User-defined name

  Surface(pugi::xml_node surf_node);

  //! Determine which side of a surface a point lies on.
  //! @param xyz[3] The 3D Cartesian coordinate of a point.
  //! @param uvw[3] A direction used to "break ties" and pick a sense when the
  //!   point is very close to the surface.
  //! @return true if the point is on the "positive" side of the surface and
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

  //! Write all information needed to reconstruct the surface to an HDF5 group.
  //! @param group_id An HDF5 group id.
  void to_hdf5(hid_t group_id) const;

protected:
  virtual void to_hdf5_inner(hid_t group_id) const = 0;
};

//==============================================================================
//! A `Surface` that supports periodic boundary conditions.
//!
//! Translational periodicity is supported for the `XPlane`, `YPlane`, `ZPlane`,
//! and `Plane` types.  Rotational periodicity is supported for
//! `XPlane`-`YPlane` pairs.
//==============================================================================

class PeriodicSurface : public Surface
{
public:
  int i_periodic{C_NONE};    //!< Index of corresponding periodic surface

  PeriodicSurface(pugi::xml_node surf_node);

  //! Translate a particle onto this surface from a periodic partner surface.
  //! @param other A pointer to the partner surface in this periodic BC.
  //! @param xyz[3] A point on the partner surface that will be translated onto
  //!   this surface.
  //! @param uvw[3] A direction that will be rotated for systems with rotational
  //!   periodicity.
  //! @return true if this surface and its partner make a rotationally-periodic
  //!   boundary condition.
  virtual bool periodic_translate(PeriodicSurface *other, double xyz[3],
                                  double uvw[3]) const = 0;

  //! Get the bounding box for this surface.
  virtual struct BoundingBox bounding_box() const = 0;
};

//==============================================================================
//! A plane perpendicular to the x-axis.
//
//! The plane is described by the equation \f$x - x_0 = 0\f$
//==============================================================================

class SurfaceXPlane : public PeriodicSurface
{
  double x0;
public:
  SurfaceXPlane(pugi::xml_node surf_node);
  double evaluate(const double xyz[3]) const;
  double distance(const double xyz[3], const double uvw[3], bool coincident)
         const;
  void normal(const double xyz[3], double uvw[3]) const;
  void to_hdf5_inner(hid_t group_id) const;
  bool periodic_translate(PeriodicSurface *other, double xyz[3], double uvw[3])
       const;
  struct BoundingBox bounding_box() const;
};

//==============================================================================
//! A plane perpendicular to the y-axis.
//
//! The plane is described by the equation \f$y - y_0 = 0\f$
//==============================================================================

class SurfaceYPlane : public PeriodicSurface
{
  double y0;
public:
  SurfaceYPlane(pugi::xml_node surf_node);
  double evaluate(const double xyz[3]) const;
  double distance(const double xyz[3], const double uvw[3],
                  bool coincident) const;
  void normal(const double xyz[3], double uvw[3]) const;
  void to_hdf5_inner(hid_t group_id) const;
  bool periodic_translate(PeriodicSurface *other, double xyz[3], double uvw[3])
       const;
  struct BoundingBox bounding_box() const;
};

//==============================================================================
//! A plane perpendicular to the z-axis.
//
//! The plane is described by the equation \f$z - z_0 = 0\f$
//==============================================================================

class SurfaceZPlane : public PeriodicSurface
{
  double z0;
public:
  SurfaceZPlane(pugi::xml_node surf_node);
  double evaluate(const double xyz[3]) const;
  double distance(const double xyz[3], const double uvw[3],
                  bool coincident) const;
  void normal(const double xyz[3], double uvw[3]) const;
  void to_hdf5_inner(hid_t group_id) const;
  bool periodic_translate(PeriodicSurface *other, double xyz[3], double uvw[3])
       const;
  struct BoundingBox bounding_box() const;
};

//==============================================================================
//! A general plane.
//
//! The plane is described by the equation \f$A x + B y + C z - D = 0\f$
//==============================================================================

class SurfacePlane : public PeriodicSurface
{
  double A, B, C, D;
public:
  SurfacePlane(pugi::xml_node surf_node);
  double evaluate(const double xyz[3]) const;
  double distance(const double xyz[3], const double uvw[3], bool coincident)
         const;
  void normal(const double xyz[3], double uvw[3]) const;
  void to_hdf5_inner(hid_t group_id) const;
  bool periodic_translate(PeriodicSurface *other, double xyz[3], double uvw[3])
       const;
  struct BoundingBox bounding_box() const;
};

//==============================================================================
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
  void to_hdf5_inner(hid_t group_id) const;
};

//==============================================================================
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
  void to_hdf5_inner(hid_t group_id) const;
};

//==============================================================================
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
  void to_hdf5_inner(hid_t group_id) const;
};

//==============================================================================
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
  void to_hdf5_inner(hid_t group_id) const;
};

//==============================================================================
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
  void to_hdf5_inner(hid_t group_id) const;
};

//==============================================================================
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
  void to_hdf5_inner(hid_t group_id) const;
};

//==============================================================================
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
  void to_hdf5_inner(hid_t group_id) const;
};

//==============================================================================
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
  void to_hdf5_inner(hid_t group_id) const;
};

//==============================================================================
// Fortran compatibility functions
//==============================================================================

extern "C" bool
surface_sense(int surf_ind, double xyz[3], double uvw[3])
{
  return surfaces_c[surf_ind]->sense(xyz, uvw);
}

extern "C" void
surface_reflect(int surf_ind, double xyz[3], double uvw[3])
{
  surfaces_c[surf_ind]->reflect(xyz, uvw);
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

extern "C" void
surface_to_hdf5(int surf_ind, hid_t group)
{
  surfaces_c[surf_ind]->to_hdf5(group);
}

extern "C" bool
surface_periodic(int surf_ind1, double xyz[3], double uvw[3])
{
  // Hopefully this function has only been called for a pair of surfaces that
  // support periodic BCs (checking should have been done when reading the
  // geometry XML).  Downcast the surfaces to the PeriodicSurface type so we
  // can call the periodic_translate method.
  Surface *surf1_gen = surfaces_c[surf_ind1];
  PeriodicSurface *surf1 = dynamic_cast<PeriodicSurface *>(surf1_gen);
  Surface *surf2_gen = surfaces_c[surf1->i_periodic];
  PeriodicSurface *surf2 = dynamic_cast<PeriodicSurface *>(surf2_gen);

  // Call the type-bound methods.
  return surf2->periodic_translate(surf1, xyz, uvw);
}

extern "C" int
surface_i_periodic(int surf_ind)
{
  PeriodicSurface *surf = dynamic_cast<PeriodicSurface *>(surfaces_c[surf_ind]);
  return surf->i_periodic;
}

#endif // SURFACE_H
