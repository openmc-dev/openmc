#ifndef OPENMC_SURFACE_H
#define OPENMC_SURFACE_H

#include <limits> // For numeric_limits
#include <list>
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

#include "tpms.h"
#include "tpms_functions.h"

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
  bool surf_source_ {false}; //!< Activate source banking for the surface?

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

  //! Determine if the surface has a special type and needs specific ray-tracing
  //! treatment. \return the name of the special treatment.
  virtual bool is_tpms() const { return false; };

  //! Determine the direction of a ray reflected from the surface.
  //! \param[in] r The point at which the ray is incident.
  //! \param[in] u Incident direction of the ray
  //! \param[inout] p Pointer to the particle. Only DAGMC uses this.
  //! \return Outgoing direction of the ray
  virtual Direction reflect(
    Position r, Direction u, GeometryState* p = nullptr) const;

  virtual Direction diffuse_reflect(
    Position r, Direction u, uint64_t* seed) const;

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
  virtual double distance(Position r, Direction u, bool coincident) const
  {
    return 0.;
  };
  virtual double distance(
    Position r, Direction u, bool coincident, double max_range) const
  {
    return 0.;
  };

  //! Compute the local outward normal direction of the surface.
  //! \param r A 3D Cartesian coordinate.
  //! \return Normal direction
  virtual Direction normal(Position r) const = 0;

  //! Write all information needed to reconstruct the surface to an HDF5 group.
  //! \param group_id An HDF5 group id.
  void to_hdf5(hid_t group_id) const;

  //! Get the BoundingBox for this surface.
  virtual BoundingBox bounding_box(bool /*pos_side*/) const { return {}; }

  /* Must specify if this is a CSG or DAGMC-type surface. Only
   * the DAGMC surface should return the DAG type geometry, so
   * by default, this returns the CSG. The main difference is that
   * if the geom_type is found to be DAG in the geometry handling code,
   * some DAGMC-specific operations get carried out like resetting
   * the particle's intersection history when necessary.
   */
  virtual GeometryType geom_type() const { return GeometryType::CSG; }

protected:
  virtual void to_hdf5_inner(hid_t group_id) const = 0;
};

//==============================================================================
//! A plane perpendicular to the x-axis.
//
//! The plane is described by the equation \f$x - x_0 = 0\f$
//==============================================================================

class SurfaceXPlane : public Surface {
public:
  explicit SurfaceXPlane(pugi::xml_node surf_node);
  double evaluate(Position r) const override;
  double distance(Position r, Direction u, bool coincident) const override;
  Direction normal(Position r) const override;
  void to_hdf5_inner(hid_t group_id) const override;
  BoundingBox bounding_box(bool pos_side) const override;

  double x0_;
};

//==============================================================================
//! A plane perpendicular to the y-axis.
//
//! The plane is described by the equation \f$y - y_0 = 0\f$
//==============================================================================

class SurfaceYPlane : public Surface {
public:
  explicit SurfaceYPlane(pugi::xml_node surf_node);
  double evaluate(Position r) const override;
  double distance(Position r, Direction u, bool coincident) const override;
  Direction normal(Position r) const override;
  void to_hdf5_inner(hid_t group_id) const override;
  BoundingBox bounding_box(bool pos_side) const override;

  double y0_;
};

//==============================================================================
//! A plane perpendicular to the z-axis.
//
//! The plane is described by the equation \f$z - z_0 = 0\f$
//==============================================================================

class SurfaceZPlane : public Surface {
public:
  explicit SurfaceZPlane(pugi::xml_node surf_node);
  double evaluate(Position r) const override;
  double distance(Position r, Direction u, bool coincident) const override;
  Direction normal(Position r) const override;
  void to_hdf5_inner(hid_t group_id) const override;
  BoundingBox bounding_box(bool pos_side) const override;

  double z0_;
};

//==============================================================================
//! A general plane.
//
//! The plane is described by the equation \f$A x + B y + C z - D = 0\f$
//==============================================================================

class SurfacePlane : public Surface {
public:
  explicit SurfacePlane(pugi::xml_node surf_node);
  double evaluate(Position r) const override;
  double distance(Position r, Direction u, bool coincident) const override;
  Direction normal(Position r) const override;
  void to_hdf5_inner(hid_t group_id) const override;

  double A_, B_, C_, D_;
};

//==============================================================================
//! A cylinder aligned along the x-axis.
//
//! The cylinder is described by the equation
//! \f$(y - y_0)^2 + (z - z_0)^2 - R^2 = 0\f$
//==============================================================================

class SurfaceXCylinder : public Surface {
public:
  explicit SurfaceXCylinder(pugi::xml_node surf_node);
  double evaluate(Position r) const override;
  double distance(Position r, Direction u, bool coincident) const override;
  Direction normal(Position r) const override;
  void to_hdf5_inner(hid_t group_id) const override;
  BoundingBox bounding_box(bool pos_side) const override;

  double y0_, z0_, radius_;
};

//==============================================================================
//! A cylinder aligned along the y-axis.
//
//! The cylinder is described by the equation
//! \f$(x - x_0)^2 + (z - z_0)^2 - R^2 = 0\f$
//==============================================================================

class SurfaceYCylinder : public Surface {
public:
  explicit SurfaceYCylinder(pugi::xml_node surf_node);
  double evaluate(Position r) const override;
  double distance(Position r, Direction u, bool coincident) const override;
  Direction normal(Position r) const override;
  void to_hdf5_inner(hid_t group_id) const override;
  BoundingBox bounding_box(bool pos_side) const override;

  double x0_, z0_, radius_;
};

//==============================================================================
//! A cylinder aligned along the z-axis.
//
//! The cylinder is described by the equation
//! \f$(x - x_0)^2 + (y - y_0)^2 - R^2 = 0\f$
//==============================================================================

class SurfaceZCylinder : public Surface {
public:
  explicit SurfaceZCylinder(pugi::xml_node surf_node);
  double evaluate(Position r) const override;
  double distance(Position r, Direction u, bool coincident) const override;
  Direction normal(Position r) const override;
  void to_hdf5_inner(hid_t group_id) const override;
  BoundingBox bounding_box(bool pos_side) const override;

  double x0_, y0_, radius_;
};

//==============================================================================
//! A sphere.
//
//! The cylinder is described by the equation
//! \f$(x - x_0)^2 + (y - y_0)^2 + (z - z_0)^2 - R^2 = 0\f$
//==============================================================================

class SurfaceSphere : public Surface {
public:
  explicit SurfaceSphere(pugi::xml_node surf_node);
  double evaluate(Position r) const override;
  double distance(Position r, Direction u, bool coincident) const override;
  Direction normal(Position r) const override;
  void to_hdf5_inner(hid_t group_id) const override;
  BoundingBox bounding_box(bool pos_side) const override;

  double x0_, y0_, z0_, radius_;
};

//==============================================================================
//! A cone aligned along the x-axis.
//
//! The cylinder is described by the equation
//! \f$(y - y_0)^2 + (z - z_0)^2 - R^2 (x - x_0)^2 = 0\f$
//==============================================================================

class SurfaceXCone : public Surface {
public:
  explicit SurfaceXCone(pugi::xml_node surf_node);
  double evaluate(Position r) const override;
  double distance(Position r, Direction u, bool coincident) const override;
  Direction normal(Position r) const override;
  void to_hdf5_inner(hid_t group_id) const override;

  double x0_, y0_, z0_, radius_sq_;
};

//==============================================================================
//! A cone aligned along the y-axis.
//
//! The cylinder is described by the equation
//! \f$(x - x_0)^2 + (z - z_0)^2 - R^2 (y - y_0)^2 = 0\f$
//==============================================================================

class SurfaceYCone : public Surface {
public:
  explicit SurfaceYCone(pugi::xml_node surf_node);
  double evaluate(Position r) const override;
  double distance(Position r, Direction u, bool coincident) const override;
  Direction normal(Position r) const override;
  void to_hdf5_inner(hid_t group_id) const override;

  double x0_, y0_, z0_, radius_sq_;
};

//==============================================================================
//! A cone aligned along the z-axis.
//
//! The cylinder is described by the equation
//! \f$(x - x_0)^2 + (y - y_0)^2 - R^2 (z - z_0)^2 = 0\f$
//==============================================================================

class SurfaceZCone : public Surface {
public:
  explicit SurfaceZCone(pugi::xml_node surf_node);
  double evaluate(Position r) const override;
  double distance(Position r, Direction u, bool coincident) const override;
  Direction normal(Position r) const override;
  void to_hdf5_inner(hid_t group_id) const override;

  double x0_, y0_, z0_, radius_sq_;
};

//==============================================================================
//! A general surface described by a quadratic equation.
//
//! \f$A x^2 + B y^2 + C z^2 + D x y + E y z + F x z + G x + H y + J z + K =
//! 0\f$
//==============================================================================

class SurfaceQuadric : public Surface {
public:
  explicit SurfaceQuadric(pugi::xml_node surf_node);
  double evaluate(Position r) const override;
  double distance(Position r, Direction u, bool coincident) const override;
  Direction normal(Position r) const override;
  void to_hdf5_inner(hid_t group_id) const override;

  // Ax^2 + By^2 + Cz^2 + Dxy + Eyz + Fxz + Gx + Hy + Jz + K = 0
  double A_, B_, C_, D_, E_, F_, G_, H_, J_, K_;
};

//==============================================================================
//! Triply Periodic Minimal Surfaces (TPMS)
//
//! Surfaces defined by an implicit function f(x,y,z) = c with c, an isovalue.
//! Various equations exists :
//! * Schwarz-P : cos(l*x) + cos(l*y) + cos(l*z) = c
//! * Gyroid : sin(l*x)cos(l*z) + sin(l*y)cos(l*x) + sin(l*z)cos(l*y) = c
//! * Diamond : sin(l*x)sin(l*y)sin(l*z) + sin(l*x)cos(l*y)cos(l*z) +
//!             cos(l*x)sin(l*y)cos(l*z) + cos(l*x)cos(l*y)sin(l*z) = c
//! with l=2pi/L, L is defined as the pitch of the TPMS lattice.
//! L and c are constants for this class.
//==============================================================================

class SurfaceTPMS : public CSGSurface {
public:
  explicit SurfaceTPMS(pugi::xml_node surf_node);
  bool is_tpms() const { return true; };
  double evaluate(Position r) const;
  double distance(
    Position r, Direction u, bool coincident, double max_range) const;
  Direction normal(Position r) const;
  void to_hdf5_inner(hid_t group_id) const;

  double isovalue, pitch;
  double x0, y0, z0;
  double a, b, c, d, e, f, g, h, i;
  std::string surface_type;

private:
  std::array<std::string, 3> tpms_types = {"Schwarz_P", "Gyroid", "Diamond"};
};

//==============================================================================
//! Interpolated Triply Periodic Minimal Surfaces (TPMS)
//
//! Same surfaces as above but with isovalue c and pitch L defined as a function
//! of x,y and z.
//==============================================================================

class SurfaceFunctionTPMS : public CSGSurface {
public:
  explicit SurfaceFunctionTPMS(pugi::xml_node surf_node);
  ~SurfaceFunctionTPMS();
  bool is_tpms() const { return true; };
  double evaluate(Position r) const;
  double distance(
    Position r, Direction u, bool coincident, double max_range) const;
  Direction normal(Position r) const;
  void to_hdf5_inner(hid_t group_id) const;

  unique_ptr<FunctionForTPMS> fPitch;
  unique_ptr<FunctionForTPMS> fIsovalue;
  double x0, y0, z0;
  double a, b, c, d, e, f, g, h, i;
  std::string surface_type;
  std::string function_type;

private:
  std::array<std::string, 3> tpms_types = {"Schwarz_P", "Gyroid", "Diamond"};
};

//==============================================================================
//! A toroidal surface described by the quartic torus lies in the x direction
//
//! \f$(x-x_0)^2/B^2 + (\sqrt{(y-y_0)^2 + (z-z_0)^2} - A)^2/C^2 -1 \f$
//==============================================================================

class SurfaceXTorus : public Surface {
public:
  explicit SurfaceXTorus(pugi::xml_node surf_node);
  double evaluate(Position r) const override;
  double distance(Position r, Direction u, bool coincident) const override;
  Direction normal(Position r) const override;
  void to_hdf5_inner(hid_t group_id) const override;

  double x0_, y0_, z0_, A_, B_, C_;
};

//==============================================================================
//! A toroidal surface described by the quartic torus lies in the y direction
//
//! \f$(y-y_0)^2/B^2 + (\sqrt{(x-x_0)^2 + (z-z_0)^2} - A)^2/C^2 -1 \f$
//==============================================================================

class SurfaceYTorus : public Surface {
public:
  explicit SurfaceYTorus(pugi::xml_node surf_node);
  double evaluate(Position r) const override;
  double distance(Position r, Direction u, bool coincident) const override;
  Direction normal(Position r) const override;
  void to_hdf5_inner(hid_t group_id) const override;

  double x0_, y0_, z0_, A_, B_, C_;
};

//==============================================================================
//! A toroidal surface described by the quartic torus lies in the z direction
//
//! \f$(z-z_0)^2/B^2 + (\sqrt{(x-x_0)^2 + (y-y_0)^2} - A)^2/C^2 -1 \f$
//==============================================================================

class SurfaceZTorus : public Surface {
public:
  explicit SurfaceZTorus(pugi::xml_node surf_node);
  double evaluate(Position r) const override;
  double distance(Position r, Direction u, bool coincident) const override;
  Direction normal(Position r) const override;
  void to_hdf5_inner(hid_t group_id) const override;

  double x0_, y0_, z0_, A_, B_, C_;
};

//==============================================================================
// Non-member functions
//==============================================================================

void read_surfaces(pugi::xml_node node);

void free_memory_surfaces();

} // namespace openmc
#endif // OPENMC_SURFACE_H
