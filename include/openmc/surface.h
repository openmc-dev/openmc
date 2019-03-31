#ifndef OPENMC_SURFACE_H
#define OPENMC_SURFACE_H

#include <memory>  // for unique_ptr
#include <limits>  // For numeric_limits
#include <string>
#include <unordered_map>
#include <vector>

#include "hdf5.h"
#include "pugixml.hpp"

#include "openmc/constants.h"
#include "openmc/position.h"

#ifdef DAGMC
#include "DagMC.hpp"
#endif

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

class Surface;

namespace model {
  extern std::vector<std::unique_ptr<Surface>> surfaces;
  extern std::unordered_map<int, int> surface_map;
} // namespace model

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
  int id_;                    //!< Unique ID
  int bc_;                    //!< Boundary condition
  std::string name_;          //!< User-defined name

  explicit Surface(pugi::xml_node surf_node);
  Surface();

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
  virtual void to_hdf5(hid_t group_id) const = 0;

};

class CSGSurface : public Surface
{
public:
  explicit CSGSurface(pugi::xml_node surf_node);
  CSGSurface();

  void to_hdf5(hid_t group_id) const;

protected:
  virtual void to_hdf5_inner(hid_t group_id) const = 0;
};

//==============================================================================
//! A `Surface` representing a DAGMC-based surface in DAGMC.
//==============================================================================
#ifdef DAGMC
class DAGSurface : public Surface
{
public:
  moab::DagMC* dagmc_ptr_;
  DAGSurface();
  int32_t dag_index_;

  double evaluate(Position r) const;
  double distance(Position r, Direction u, bool coincident) const;
  Direction normal(Position r) const;
  //! Get the bounding box of this surface.
  BoundingBox bounding_box() const;

  void to_hdf5(hid_t group_id) const;
};
#endif
//==============================================================================
//! A `Surface` that supports periodic boundary conditions.
//!
//! Translational periodicity is supported for the `XPlane`, `YPlane`, `ZPlane`,
//! and `Plane` types.  Rotational periodicity is supported for
//! `XPlane`-`YPlane` pairs.
//==============================================================================

class PeriodicSurface : public CSGSurface
{
public:
  int i_periodic_{C_NONE};    //!< Index of corresponding periodic surface

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
public:
  explicit SurfaceXPlane(pugi::xml_node surf_node);
  double evaluate(Position r) const;
  double distance(Position r, Direction u, bool coincident) const;
  Direction normal(Position r) const;
  void to_hdf5_inner(hid_t group_id) const;
  bool periodic_translate(const PeriodicSurface* other, Position& r,
                          Direction& u) const;
  BoundingBox bounding_box() const;

  double x0_;
};

//==============================================================================
//! A plane perpendicular to the y-axis.
//
//! The plane is described by the equation \f$y - y_0 = 0\f$
//==============================================================================

class SurfaceYPlane : public PeriodicSurface
{
public:
  explicit SurfaceYPlane(pugi::xml_node surf_node);
  double evaluate(Position r) const;
  double distance(Position r, Direction u, bool coincident) const;
  Direction normal(Position r) const;
  void to_hdf5_inner(hid_t group_id) const;
  bool periodic_translate(const PeriodicSurface* other, Position& r,
                          Direction& u) const;
  BoundingBox bounding_box() const;

  double y0_;
};

//==============================================================================
//! A plane perpendicular to the z-axis.
//
//! The plane is described by the equation \f$z - z_0 = 0\f$
//==============================================================================

class SurfaceZPlane : public PeriodicSurface
{
public:
  explicit SurfaceZPlane(pugi::xml_node surf_node);
  double evaluate(Position r) const;
  double distance(Position r, Direction u, bool coincident) const;
  Direction normal(Position r) const;
  void to_hdf5_inner(hid_t group_id) const;
  bool periodic_translate(const PeriodicSurface* other, Position& r,
                          Direction& u) const;
  BoundingBox bounding_box() const;

  double z0_;
};

//==============================================================================
//! A general plane.
//
//! The plane is described by the equation \f$A x + B y + C z - D = 0\f$
//==============================================================================

class SurfacePlane : public PeriodicSurface
{
public:
  explicit SurfacePlane(pugi::xml_node surf_node);
  double evaluate(Position r) const;
  double distance(Position r, Direction u, bool coincident) const;
  Direction normal(Position r) const;
  void to_hdf5_inner(hid_t group_id) const;
  bool periodic_translate(const PeriodicSurface* other, Position& r,
                          Direction& u) const;
  BoundingBox bounding_box() const;

  double A_, B_, C_, D_;
};

//==============================================================================
//! A cylinder aligned along the x-axis.
//
//! The cylinder is described by the equation
//! \f$(y - y_0)^2 + (z - z_0)^2 - R^2 = 0\f$
//==============================================================================

class SurfaceXCylinder : public CSGSurface
{
public:
  explicit SurfaceXCylinder(pugi::xml_node surf_node);
  double evaluate(Position r) const;
  double distance(Position r, Direction u, bool coincident) const;
  Direction normal(Position r) const;
  void to_hdf5_inner(hid_t group_id) const;

  double y0_, z0_, radius_;
};

//==============================================================================
//! A cylinder aligned along the y-axis.
//
//! The cylinder is described by the equation
//! \f$(x - x_0)^2 + (z - z_0)^2 - R^2 = 0\f$
//==============================================================================

class SurfaceYCylinder : public CSGSurface
{
public:
  explicit SurfaceYCylinder(pugi::xml_node surf_node);
  double evaluate(Position r) const;
  double distance(Position r, Direction u, bool coincident) const;
  Direction normal(Position r) const;
  void to_hdf5_inner(hid_t group_id) const;

  double x0_, z0_, radius_;
};

//==============================================================================
//! A cylinder aligned along the z-axis.
//
//! The cylinder is described by the equation
//! \f$(x - x_0)^2 + (y - y_0)^2 - R^2 = 0\f$
//==============================================================================

class SurfaceZCylinder : public CSGSurface
{
public:
  explicit SurfaceZCylinder(pugi::xml_node surf_node);
  double evaluate(Position r) const;
  double distance(Position r, Direction u, bool coincident) const;
  Direction normal(Position r) const;
  void to_hdf5_inner(hid_t group_id) const;

  double x0_, y0_, radius_;
};

//==============================================================================
//! A sphere.
//
//! The cylinder is described by the equation
//! \f$(x - x_0)^2 + (y - y_0)^2 + (z - z_0)^2 - R^2 = 0\f$
//==============================================================================

class SurfaceSphere : public CSGSurface
{
public:
  explicit SurfaceSphere(pugi::xml_node surf_node);
  double evaluate(Position r) const;
  double distance(Position r, Direction u, bool coincident) const;
  Direction normal(Position r) const;
  void to_hdf5_inner(hid_t group_id) const;

  double x0_, y0_, z0_, radius_;
};

//==============================================================================
//! A cone aligned along the x-axis.
//
//! The cylinder is described by the equation
//! \f$(y - y_0)^2 + (z - z_0)^2 - R^2 (x - x_0)^2 = 0\f$
//==============================================================================

class SurfaceXCone : public CSGSurface
{
public:
  explicit SurfaceXCone(pugi::xml_node surf_node);
  double evaluate(Position r) const;
  double distance(Position r, Direction u, bool coincident) const;
  Direction normal(Position r) const;
  void to_hdf5_inner(hid_t group_id) const;

  double x0_, y0_, z0_, radius_sq_;
};

//==============================================================================
//! A cone aligned along the y-axis.
//
//! The cylinder is described by the equation
//! \f$(x - x_0)^2 + (z - z_0)^2 - R^2 (y - y_0)^2 = 0\f$
//==============================================================================

class SurfaceYCone : public CSGSurface
{
public:
  explicit SurfaceYCone(pugi::xml_node surf_node);
  double evaluate(Position r) const;
  double distance(Position r, Direction u, bool coincident) const;
  Direction normal(Position r) const;
  void to_hdf5_inner(hid_t group_id) const;

  double x0_, y0_, z0_, radius_sq_;
};

//==============================================================================
//! A cone aligned along the z-axis.
//
//! The cylinder is described by the equation
//! \f$(x - x_0)^2 + (y - y_0)^2 - R^2 (z - z_0)^2 = 0\f$
//==============================================================================

class SurfaceZCone : public CSGSurface
{
public:
  explicit SurfaceZCone(pugi::xml_node surf_node);
  double evaluate(Position r) const;
  double distance(Position r, Direction u, bool coincident) const;
  Direction normal(Position r) const;
  void to_hdf5_inner(hid_t group_id) const;

  double x0_, y0_, z0_, radius_sq_;
};

//==============================================================================
//! A general surface described by a quadratic equation.
//
//! \f$A x^2 + B y^2 + C z^2 + D x y + E y z + F x z + G x + H y + J z + K = 0\f$
//==============================================================================

class SurfaceQuadric : public CSGSurface
{
public:
  explicit SurfaceQuadric(pugi::xml_node surf_node);
  double evaluate(Position r) const;
  double distance(Position r, Direction u, bool coincident) const;
  Direction normal(Position r) const;
  void to_hdf5_inner(hid_t group_id) const;

  // Ax^2 + By^2 + Cz^2 + Dxy + Eyz + Fxz + Gx + Hy + Jz + K = 0
  double A_, B_, C_, D_, E_, F_, G_, H_, J_, K_;
};

//==============================================================================
// Non-member functions
//==============================================================================

void read_surfaces(pugi::xml_node node);

void free_memory_surfaces();

} // namespace openmc
#endif // OPENMC_SURFACE_H
