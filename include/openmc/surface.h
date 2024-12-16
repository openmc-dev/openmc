#ifndef OPENMC_SURFACE_H
#define OPENMC_SURFACE_H

#include <limits> // For numeric_limits
#include <string>
#include <unordered_map>

#include "hdf5.h"
#include "pugixml.hpp"

#include "openmc/boundary_condition.h"
#include "openmc/bounding_box.h"
#include "openmc/constants.h"
#include "openmc/memory.h" // for unique_ptr
#include "openmc/particle.h"
#include "openmc/position.h"
#include "openmc/vector.h"

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

class Surface;

namespace model {
extern std::unordered_map<int, int> surface_map;
extern vector<unique_ptr<Surface>> surfaces;
} // namespace model

//==============================================================================
//! A geometry primitive used to define regions of 3D space.
//==============================================================================

class Surface {
public:
  int id_;                           //!< Unique ID
  std::string name_;                 //!< User-defined name
  unique_ptr<BoundaryCondition> bc_; //!< Boundary condition
  GeometryType geom_type_;           //!< Geometry type indicator (CSG or DAGMC)
  bool surf_source_ {false}; //!< Activate source banking for the surface?

  explicit Surface(pugi::xml_node surf_node);
  Surface();

  virtual ~Surface() {}

  //! Determine which side of a surface a point lies on.
  //! \param r The 3D Cartesian coordinate of a point.
  //! \param u A direction used to "break ties" and pick a sense when the
  //!   point is very close to the surface.
  //! \param t The time for the evaluation
  //! \param speed The point speed to "break ties" with moving surface.
  //! \return true if the point is on the "positive" side of the surface and
  //!   false otherwise.
  bool sense(Position r, Direction u, double t, double speed) const;

  //! Determine the direction of a ray reflected from the surface.
  //! \param[in] r The point at which the ray is incident.
  //! \param[in] u Incident direction of the ray
  //! \param[inout] p Pointer to the particle. Only DAGMC uses this.
  //! \return Outgoing direction of the ray
  virtual Direction reflect(
    Position r, Direction u, GeometryState* p = nullptr) const;

  virtual Direction diffuse_reflect(
    Position r, Direction u, uint64_t* seed, GeometryState* p = nullptr) const;

  //! Evaluate the equation describing the surface.
  //!
  //! Surfaces can be described by some function f(x, y, z) = 0. This member
  //! function evaluates that mathematical function. The time variable t is
  //! needed for evaluation of moving surfaces.
  //! \param r A 3D Cartesian coordinate.
  //! \param t The time for the evaluation.
  virtual double evaluate(Position r, double t) const = 0;

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

  //! Compute the dot product of the local outward normal direction of the
  //! surface to a geometry coordinate.
  //! \param r A 3D Cartesian coordinate.
  //! \param u The direction of the ray.
  //! \param t The time for the evaluation.
  //! \param speed The speed of the particle.
  //! \return The dot product
  virtual double dot_normal(
    Position r, Direction u, double t, double speed) const = 0;

  //! Write all information needed to reconstruct the surface to an HDF5 group.
  //! \param group_id An HDF5 group id.
  void to_hdf5(hid_t group_id) const;

  //! Get the BoundingBox for this surface.
  virtual BoundingBox bounding_box(bool /*pos_side*/) const { return {}; }

protected:
  virtual void to_hdf5_inner(hid_t group_id) const = 0;
};

class CSGSurface : public Surface {
public:
  //!< Moving surface parameters
  bool moving_ {false};
  vector<double> moving_time_grid_;      // Time grid points [0.0, ..., INFTY]
  vector<Position> moving_translations_; // Translations at the time grid points
  vector<Position> moving_velocities_;   // Velocities within time bins

  explicit CSGSurface(pugi::xml_node surf_node);
  CSGSurface();
  double evaluate(Position r, double t) const override;
  double dot_normal(
    Position r, Direction u, double t, double speed) const override;

protected:
  //! Static CSG surface evaluation.
  virtual double _evaluate(Position r) const = 0;
};

//==============================================================================
//! A plane perpendicular to the x-axis.
//
//! The plane is described by the equation \f$x - x_0 = 0\f$
//==============================================================================

class SurfaceXPlane : public CSGSurface {
public:
  explicit SurfaceXPlane(pugi::xml_node surf_node);
  double distance(Position r, Direction u, bool coincident) const override;
  Direction normal(Position r) const override;
  void to_hdf5_inner(hid_t group_id) const override;
  BoundingBox bounding_box(bool pos_side) const override;

  double x0_;

protected:
  double _evaluate(Position r) const override;
};

//==============================================================================
//! A plane perpendicular to the y-axis.
//
//! The plane is described by the equation \f$y - y_0 = 0\f$
//==============================================================================

class SurfaceYPlane : public CSGSurface {
public:
  explicit SurfaceYPlane(pugi::xml_node surf_node);
  double distance(Position r, Direction u, bool coincident) const override;
  Direction normal(Position r) const override;
  void to_hdf5_inner(hid_t group_id) const override;
  BoundingBox bounding_box(bool pos_side) const override;

  double y0_;

protected:
  double _evaluate(Position r) const override;
};

//==============================================================================
//! A plane perpendicular to the z-axis.
//
//! The plane is described by the equation \f$z - z_0 = 0\f$
//==============================================================================

class SurfaceZPlane : public CSGSurface {
public:
  explicit SurfaceZPlane(pugi::xml_node surf_node);
  double distance(Position r, Direction u, bool coincident) const override;
  Direction normal(Position r) const override;
  void to_hdf5_inner(hid_t group_id) const override;
  BoundingBox bounding_box(bool pos_side) const override;

  double z0_;

protected:
  double _evaluate(Position r) const override;
};

//==============================================================================
//! A general plane.
//
//! The plane is described by the equation \f$A x + B y + C z - D = 0\f$
//==============================================================================

class SurfacePlane : public CSGSurface {
public:
  explicit SurfacePlane(pugi::xml_node surf_node);
  double distance(Position r, Direction u, bool coincident) const override;
  Direction normal(Position r) const override;
  void to_hdf5_inner(hid_t group_id) const override;

  double A_, B_, C_, D_;

protected:
  double _evaluate(Position r) const override;
};

//==============================================================================
//! A cylinder aligned along the x-axis.
//
//! The cylinder is described by the equation
//! \f$(y - y_0)^2 + (z - z_0)^2 - R^2 = 0\f$
//==============================================================================

class SurfaceXCylinder : public CSGSurface {
public:
  explicit SurfaceXCylinder(pugi::xml_node surf_node);
  double distance(Position r, Direction u, bool coincident) const override;
  Direction normal(Position r) const override;
  void to_hdf5_inner(hid_t group_id) const override;
  BoundingBox bounding_box(bool pos_side) const override;

  double y0_, z0_, radius_;

protected:
  double _evaluate(Position r) const override;
};

//==============================================================================
//! A cylinder aligned along the y-axis.
//
//! The cylinder is described by the equation
//! \f$(x - x_0)^2 + (z - z_0)^2 - R^2 = 0\f$
//==============================================================================

class SurfaceYCylinder : public CSGSurface {
public:
  explicit SurfaceYCylinder(pugi::xml_node surf_node);
  double distance(Position r, Direction u, bool coincident) const override;
  Direction normal(Position r) const override;
  void to_hdf5_inner(hid_t group_id) const override;
  BoundingBox bounding_box(bool pos_side) const override;

  double x0_, z0_, radius_;

protected:
  double _evaluate(Position r) const override;
};

//==============================================================================
//! A cylinder aligned along the z-axis.
//
//! The cylinder is described by the equation
//! \f$(x - x_0)^2 + (y - y_0)^2 - R^2 = 0\f$
//==============================================================================

class SurfaceZCylinder : public CSGSurface {
public:
  explicit SurfaceZCylinder(pugi::xml_node surf_node);
  double distance(Position r, Direction u, bool coincident) const override;
  Direction normal(Position r) const override;
  void to_hdf5_inner(hid_t group_id) const override;
  BoundingBox bounding_box(bool pos_side) const override;

  double x0_, y0_, radius_;

protected:
  double _evaluate(Position r) const override;
};

//==============================================================================
//! A sphere.
//
//! The cylinder is described by the equation
//! \f$(x - x_0)^2 + (y - y_0)^2 + (z - z_0)^2 - R^2 = 0\f$
//==============================================================================

class SurfaceSphere : public CSGSurface {
public:
  explicit SurfaceSphere(pugi::xml_node surf_node);
  double distance(Position r, Direction u, bool coincident) const override;
  Direction normal(Position r) const override;
  void to_hdf5_inner(hid_t group_id) const override;
  BoundingBox bounding_box(bool pos_side) const override;

  double x0_, y0_, z0_, radius_;

protected:
  double _evaluate(Position r) const override;
};

//==============================================================================
//! A cone aligned along the x-axis.
//
//! The cylinder is described by the equation
//! \f$(y - y_0)^2 + (z - z_0)^2 - R^2 (x - x_0)^2 = 0\f$
//==============================================================================

class SurfaceXCone : public CSGSurface {
public:
  explicit SurfaceXCone(pugi::xml_node surf_node);
  double distance(Position r, Direction u, bool coincident) const override;
  Direction normal(Position r) const override;
  void to_hdf5_inner(hid_t group_id) const override;

  double x0_, y0_, z0_, radius_sq_;

protected:
  double _evaluate(Position r) const override;
};

//==============================================================================
//! A cone aligned along the y-axis.
//
//! The cylinder is described by the equation
//! \f$(x - x_0)^2 + (z - z_0)^2 - R^2 (y - y_0)^2 = 0\f$
//==============================================================================

class SurfaceYCone : public CSGSurface {
public:
  explicit SurfaceYCone(pugi::xml_node surf_node);
  double distance(Position r, Direction u, bool coincident) const override;
  Direction normal(Position r) const override;
  void to_hdf5_inner(hid_t group_id) const override;

  double x0_, y0_, z0_, radius_sq_;

protected:
  double _evaluate(Position r) const override;
};

//==============================================================================
//! A cone aligned along the z-axis.
//
//! The cylinder is described by the equation
//! \f$(x - x_0)^2 + (y - y_0)^2 - R^2 (z - z_0)^2 = 0\f$
//==============================================================================

class SurfaceZCone : public CSGSurface {
public:
  explicit SurfaceZCone(pugi::xml_node surf_node);
  double distance(Position r, Direction u, bool coincident) const override;
  Direction normal(Position r) const override;
  void to_hdf5_inner(hid_t group_id) const override;

  double x0_, y0_, z0_, radius_sq_;

protected:
  double _evaluate(Position r) const override;
};

//==============================================================================
//! A general surface described by a quadratic equation.
//
//! \f$A x^2 + B y^2 + C z^2 + D x y + E y z + F x z + G x + H y + J z + K =
//! 0\f$
//==============================================================================

class SurfaceQuadric : public CSGSurface {
public:
  explicit SurfaceQuadric(pugi::xml_node surf_node);
  double distance(Position r, Direction u, bool coincident) const override;
  Direction normal(Position r) const override;
  void to_hdf5_inner(hid_t group_id) const override;

  // Ax^2 + By^2 + Cz^2 + Dxy + Eyz + Fxz + Gx + Hy + Jz + K = 0
  double A_, B_, C_, D_, E_, F_, G_, H_, J_, K_;

protected:
  double _evaluate(Position r) const override;
};

//==============================================================================
//! A toroidal surface described by the quartic torus lies in the x direction
//
//! \f$(x-x_0)^2/B^2 + (\sqrt{(y-y_0)^2 + (z-z_0)^2} - A)^2/C^2 -1 \f$
//==============================================================================

class SurfaceXTorus : public CSGSurface {
public:
  explicit SurfaceXTorus(pugi::xml_node surf_node);
  double distance(Position r, Direction u, bool coincident) const override;
  Direction normal(Position r) const override;
  void to_hdf5_inner(hid_t group_id) const override;

  double x0_, y0_, z0_, A_, B_, C_;

protected:
  double _evaluate(Position r) const override;
};

//==============================================================================
//! A toroidal surface described by the quartic torus lies in the y direction
//
//! \f$(y-y_0)^2/B^2 + (\sqrt{(x-x_0)^2 + (z-z_0)^2} - A)^2/C^2 -1 \f$
//==============================================================================

class SurfaceYTorus : public CSGSurface {
public:
  explicit SurfaceYTorus(pugi::xml_node surf_node);
  double distance(Position r, Direction u, bool coincident) const override;
  Direction normal(Position r) const override;
  void to_hdf5_inner(hid_t group_id) const override;

  double x0_, y0_, z0_, A_, B_, C_;

protected:
  double _evaluate(Position r) const override;
};

//==============================================================================
//! A toroidal surface described by the quartic torus lies in the z direction
//
//! \f$(z-z_0)^2/B^2 + (\sqrt{(x-x_0)^2 + (y-y_0)^2} - A)^2/C^2 -1 \f$
//==============================================================================

class SurfaceZTorus : public CSGSurface {
public:
  explicit SurfaceZTorus(pugi::xml_node surf_node);
  double distance(Position r, Direction u, bool coincident) const override;
  Direction normal(Position r) const override;
  void to_hdf5_inner(hid_t group_id) const override;

  double x0_, y0_, z0_, A_, B_, C_;

protected:
  double _evaluate(Position r) const override;
};

//==============================================================================
// Non-member functions
//==============================================================================

void read_surfaces(pugi::xml_node node);

void free_memory_surfaces();

} // namespace openmc
#endif // OPENMC_SURFACE_H
