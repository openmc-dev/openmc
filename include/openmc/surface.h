#ifndef OPENMC_SURFACE_H
#define OPENMC_SURFACE_H

#include <map>
#include <limits>  // For numeric_limits
#include <string>
#include <vector>

#include "hdf5.h"
#include "pugixml.hpp"

#include "openmc/constants.h"
#include "openmc/position.h"


namespace openmc {

//==============================================================================
// Module constant declarations (defined in .cpp)
//==============================================================================

// TODO: Convert to enum
extern "C" const int BC_TRANSMIT;
extern "C" const int BC_VACUUM;
extern "C" const int BC_REFLECT;
extern "C" const int BC_PERIODIC;

//==============================================================================
// Global variables
//==============================================================================

extern "C" int32_t n_surfaces;

class Surface;
extern std::vector<Surface*> global_surfaces;

extern std::map<int, int> surface_map;

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
  std::string name;          //!< User-defined name

  explicit Surface(pugi::xml_node surf_node);

  virtual ~Surface() {}

  //! Determine which side of a surface a point lies on.
  //! \param r The 3D Cartesian coordinate of a point.
  //! \param u A direction used to "break ties" and pick a sense when the
  //!   point is very close to the surface.
  //! \return true if the point is on the "positive" side of the surface and
  //!   false otherwise.
  bool sense(Position r, Direction u) const;

  //! Determine the direction of a ray reflected from the surface.
  //! \param[in] r The point at which the ray is incident.
  //! \param[in] u Incident direction of the ray
  //! \return Outgoing direction of the ray
  Direction reflect(Position r, Direction u) const;

  //! Evaluate the equation describing the surface.
  //!
  //! Surfaces can be described by some function f(x, y, z) = 0.  This member
  //! function evaluates that mathematical function.
  //! \param r A 3D Cartesian coordinate.
  virtual double evaluate(Position r) const = 0;

  //! Compute the distance between a point and the surface along a ray.
  //! \param r A 3D Cartesian coordinate.
  //! \param u The direction of the ray.
  //! \param coincident A hint to the code that the given point should lie
  //!   exactly on the surface.
  virtual double distance(Position r, Direction u, bool coincident) const = 0;

  //! Compute the local outward normal direction of the surface.
  //! \param r A 3D Cartesian coordinate.
  //! \return Normal direction
  virtual Direction normal(Position r) const = 0;

  //! Write all information needed to reconstruct the surface to an HDF5 group.
  //! \param group_id An HDF5 group id.
  //TODO: this probably needs to include i_periodic for PeriodicSurface
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

  explicit PeriodicSurface(pugi::xml_node surf_node);

  //! Translate a particle onto this surface from a periodic partner surface.
  //! \param other A pointer to the partner surface in this periodic BC.
  //! \param r A point on the partner surface that will be translated onto
  //!   this surface.
  //! \param u A direction that will be rotated for systems with rotational
  //!   periodicity.
  //! \return true if this surface and its partner make a rotationally-periodic
  //!   boundary condition.
  virtual bool periodic_translate(const PeriodicSurface* other, Position& r,
                                  Direction& u) const = 0;

  //! Get the bounding box for this surface.
  virtual BoundingBox bounding_box() const = 0;
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
  explicit SurfaceXPlane(pugi::xml_node surf_node);
  double evaluate(Position r) const;
  double distance(Position r, Direction u, bool coincident) const;
  Direction normal(Position r) const;
  void to_hdf5_inner(hid_t group_id) const;
  bool periodic_translate(const PeriodicSurface* other, Position& r,
                          Direction& u) const;
  BoundingBox bounding_box() const;
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
  explicit SurfaceYPlane(pugi::xml_node surf_node);
  double evaluate(Position r) const;
  double distance(Position r, Direction u, bool coincident) const;
  Direction normal(Position r) const;
  void to_hdf5_inner(hid_t group_id) const;
  bool periodic_translate(const PeriodicSurface* other, Position& r,
                          Direction& u) const;
  BoundingBox bounding_box() const;
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
  explicit SurfaceZPlane(pugi::xml_node surf_node);
  double evaluate(Position r) const;
  double distance(Position r, Direction u, bool coincident) const;
  Direction normal(Position r) const;
  void to_hdf5_inner(hid_t group_id) const;
  bool periodic_translate(const PeriodicSurface* other, Position& r,
                          Direction& u) const;
  BoundingBox bounding_box() const;
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
  explicit SurfacePlane(pugi::xml_node surf_node);
  double evaluate(Position r) const;
  double distance(Position r, Direction u, bool coincident) const;
  Direction normal(Position r) const;
  void to_hdf5_inner(hid_t group_id) const;
  bool periodic_translate(const PeriodicSurface* other, Position& r,
                          Direction& u) const;
  BoundingBox bounding_box() const;
};

//==============================================================================
//! A cylinder aligned along the x-axis.
//
//! The cylinder is described by the equation
//! \f$(y - y_0)^2 + (z - z_0)^2 - R^2 = 0\f$
//==============================================================================

class SurfaceXCylinder : public Surface
{
  double y0, z0, radius;
public:
  explicit SurfaceXCylinder(pugi::xml_node surf_node);
  double evaluate(Position r) const;
  double distance(Position r, Direction u, bool coincident) const;
  Direction normal(Position r) const;
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
  double x0, z0, radius;
public:
  explicit SurfaceYCylinder(pugi::xml_node surf_node);
  double evaluate(Position r) const;
  double distance(Position r, Direction u, bool coincident) const;
  Direction normal(Position r) const;
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
  double x0, y0, radius;
public:
  explicit SurfaceZCylinder(pugi::xml_node surf_node);
  double evaluate(Position r) const;
  double distance(Position r, Direction u, bool coincident) const;
  Direction normal(Position r) const;
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
  double x0, y0, z0, radius;
public:
  explicit SurfaceSphere(pugi::xml_node surf_node);
  double evaluate(Position r) const;
  double distance(Position r, Direction u, bool coincident) const;
  Direction normal(Position r) const;
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
  double x0, y0, z0, radius_sq;
public:
  explicit SurfaceXCone(pugi::xml_node surf_node);
  double evaluate(Position r) const;
  double distance(Position r, Direction u, bool coincident) const;
  Direction normal(Position r) const;
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
  double x0, y0, z0, radius_sq;
public:
  explicit SurfaceYCone(pugi::xml_node surf_node);
  double evaluate(Position r) const;
  double distance(Position r, Direction u, bool coincident) const;
  Direction normal(Position r) const;
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
  double x0, y0, z0, radius_sq;
public:
  explicit SurfaceZCone(pugi::xml_node surf_node);
  double evaluate(Position r) const;
  double distance(Position r, Direction u, bool coincident) const;
  Direction normal(Position r) const;
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
  explicit SurfaceQuadric(pugi::xml_node surf_node);
  double evaluate(Position r) const;
  double distance(Position r, Direction u, bool coincident) const;
  Direction normal(Position r) const;
  void to_hdf5_inner(hid_t group_id) const;
};

//==============================================================================
// Fortran compatibility functions
//==============================================================================

extern "C" {
  Surface* surface_pointer(int surf_ind);
  int surface_id(Surface* surf);
  int surface_bc(Surface* surf);
  bool surface_sense(Surface* surf, double xyz[3], double uvw[3]);
  void surface_reflect(Surface* surf, double xyz[3], double uvw[3]);
  double surface_distance(Surface* surf, double xyz[3], double uvw[3],
                          bool coincident);
  void surface_normal(Surface* surf, double xyz[3], double uvw[3]);
  void surface_to_hdf5(Surface* surf, hid_t group);
  int surface_i_periodic(PeriodicSurface* surf);
  bool surface_periodic(PeriodicSurface* surf, PeriodicSurface* other,
                        double xyz[3], double uvw[3]);
  void free_memory_surfaces_c();
}

} // namespace openmc
#endif // OPENMC_SURFACE_H
