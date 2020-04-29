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
#include "openmc/particle.h"
#include "openmc/position.h"
#include "dagmc.h"

namespace openmc {

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
  double xmin = -INFTY;
  double xmax = INFTY;
  double ymin = -INFTY;
  double ymax = INFTY;
  double zmin = -INFTY;
  double zmax = INFTY;


  inline BoundingBox operator &(const BoundingBox& other) {
    BoundingBox result = *this;
    return result &= other;
  }

  inline BoundingBox operator |(const BoundingBox& other) {
    BoundingBox result = *this;
    return result |= other;
  }

  // intersect operator
  inline BoundingBox& operator &=(const BoundingBox& other) {
    xmin = std::max(xmin, other.xmin);
    xmax = std::min(xmax, other.xmax);
    ymin = std::max(ymin, other.ymin);
    ymax = std::min(ymax, other.ymax);
    zmin = std::max(zmin, other.zmin);
    zmax = std::min(zmax, other.zmax);
    return *this;
  }

  // union operator
  inline BoundingBox& operator |=(const BoundingBox& other) {
    xmin = std::min(xmin, other.xmin);
    xmax = std::max(xmax, other.xmax);
    ymin = std::min(ymin, other.ymin);
    ymax = std::max(ymax, other.ymax);
    zmin = std::min(zmin, other.zmin);
    zmax = std::max(zmax, other.zmax);
    return *this;
  }

};

//==============================================================================
//! A geometry primitive used to define regions of 3D space.
//==============================================================================

class Surface
{
public:

  // Types of available boundary conditions on a surface
  enum class BoundaryType {
    TRANSMIT,
    VACUUM,
    REFLECT,
    PERIODIC,
    WHITE
  };
  
  // Types of surface
  enum class SurfaceType {
        SurfaceXPlane,
        SurfaceYPlane,
        SurfaceZPlane,
        SurfacePlane,
        SurfaceXCylinder,
        SurfaceYCylinder,
        SurfaceZCylinder,
        SurfaceSphere,
        SurfaceXCone,
        SurfaceYCone,
        SurfaceZCone,
        SurfaceQuadric
  };

  int id_;                    //!< Unique ID
  BoundaryType bc_;                    //!< Boundary condition
  SurfaceType type_;                    //!< Surface type
  std::string name_;          //!< User-defined name

  // Actual values;
  int i_periodic_{C_NONE};    //!< Index of corresponding periodic surface
  double x0_, y0_, z0_, radius_, radius_sq_;
  double A_, B_, C_, D_, E_, F_, G_, H_, J_, K_;

  explicit Surface(pugi::xml_node surf_node, SurfaceType type);
  Surface();

  ~Surface() {}

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
  //! \param[inout] p Pointer to the particle
  //! \return Outgoing direction of the ray
  Direction reflect(Position r, Direction u, Particle* p) const;

  Direction diffuse_reflect(Position r, Direction u,
    uint64_t* seed) const;

  //! Evaluate the equation describing the surface.
  //!
  //! Surfaces can be described by some function f(x, y, z) = 0.  This member
  //! function evaluates that mathematical function.
  //! \param r A 3D Cartesian coordinate.
  double evaluate(Position r) const;
 double SurfaceXPlane_evaluate(Position r) const;   
 double SurfaceYPlane_evaluate(Position r) const;   
 double SurfaceZPlane_evaluate(Position r) const;   
 double SurfacePlane_evaluate(Position r) const;    
 double SurfaceXCylinder_evaluate(Position r) const;
 double SurfaceYCylinder_evaluate(Position r) const;
 double SurfaceZCylinder_evaluate(Position r) const;
 double SurfaceSphere_evaluate(Position r) const;   
 double SurfaceXCone_evaluate(Position r) const;    
 double SurfaceYCone_evaluate(Position r) const;    
 double SurfaceZCone_evaluate(Position r) const;    
 double SurfaceQuadric_evaluate(Position r) const;  

  //! Compute the distance between a point and the surface along a ray.
  //! \param r A 3D Cartesian coordinate.
  //! \param u The direction of the ray.
  //! \param coincident A hint to the code that the given point should lie
  //!   exactly on the surface.
  double distance(Position r, Direction u, bool coincident) const;
 double SurfaceXPlane_distance(Position r, Direction u, bool coincident) const;   
 double SurfaceYPlane_distance(Position r, Direction u, bool coincident) const;   
 double SurfaceZPlane_distance(Position r, Direction u, bool coincident) const;   
 double SurfacePlane_distance(Position r, Direction u, bool coincident) const;    
 double SurfaceXCylinder_distance(Position r, Direction u, bool coincident) const;
 double SurfaceYCylinder_distance(Position r, Direction u, bool coincident) const;
 double SurfaceZCylinder_distance(Position r, Direction u, bool coincident) const;
 double SurfaceSphere_distance(Position r, Direction u, bool coincident) const;   
 double SurfaceXCone_distance(Position r, Direction u, bool coincident) const;    
 double SurfaceYCone_distance(Position r, Direction u, bool coincident) const;    
 double SurfaceZCone_distance(Position r, Direction u, bool coincident) const;    
 double SurfaceQuadric_distance(Position r, Direction u, bool coincident) const;  

  //! Compute the local outward normal direction of the surface.
  //! \param r A 3D Cartesian coordinate.
  //! \return Normal direction
  Direction normal(Position r) const;
 Direction SurfaceXPlane_normal(Position r) const;   
 Direction SurfaceYPlane_normal(Position r) const;   
 Direction SurfaceZPlane_normal(Position r) const;   
 Direction SurfacePlane_normal(Position r) const;    
 Direction SurfaceXCylinder_normal(Position r) const;
 Direction SurfaceYCylinder_normal(Position r) const;
 Direction SurfaceZCylinder_normal(Position r) const;
 Direction SurfaceSphere_normal(Position r) const;   
 Direction SurfaceXCone_normal(Position r) const;    
 Direction SurfaceYCone_normal(Position r) const;    
 Direction SurfaceZCone_normal(Position r) const;    
 Direction SurfaceQuadric_normal(Position r) const;  

  //! Write all information needed to reconstruct the surface to an HDF5 group.
  //! \param group_id An HDF5 group id.
  //TODO: this probably needs to include i_periodic for PeriodicSurface
  void to_hdf5(hid_t group_id) const;

  //! Get the BoundingBox for this surface.
  // TODO: Should this be prototype only?
  //BoundingBox bounding_box(bool pos_side) const { return {}; }
  BoundingBox bounding_box(bool pos_side) const;
 BoundingBox SurfaceXPlane_bounding_box(bool pos_side) const;   
 BoundingBox SurfaceYPlane_bounding_box(bool pos_side) const;   
 BoundingBox SurfaceZPlane_bounding_box(bool pos_side) const;   
 //BoundingBox SurfacePlane_bounding_box(bool pos_side) const;    
 BoundingBox SurfaceXCylinder_bounding_box(bool pos_side) const;
 BoundingBox SurfaceYCylinder_bounding_box(bool pos_side) const;
 BoundingBox SurfaceZCylinder_bounding_box(bool pos_side) const;
 BoundingBox SurfaceSphere_bounding_box(bool pos_side) const;   
 //BouningBox SurfaceXCone_bounding_box(bool pos_side) const;    
 //BouningBox SurfaceYCone_bounding_box(bool pos_side) const;    
 //BouningBox SurfaceZCone_bounding_box(bool pos_side) const;    
 //BouningBox SurfaceQuadric_bounding_box(bool pos_side) const;  
  
  bool periodic_translate(const Surface* other, Position& r,
                                  Direction& u) const;
 bool SurfaceXPlane_periodic_translate(const Surface* other, Position& r, Direction& u) const;   
 bool SurfaceYPlane_periodic_translate(const Surface* other, Position& r, Direction& u) const;   
 bool SurfaceZPlane_periodic_translate(const Surface* other, Position& r, Direction& u) const;   
 bool SurfacePlane_periodic_translate(const Surface* other, Position& r, Direction& u) const;    
  
protected:
  void to_hdf5_inner(hid_t group_id) const;
 void SurfaceXPlane_to_hdf5_inner(hid_t group_id) const;   
 void SurfaceYPlane_to_hdf5_inner(hid_t group_id) const;   
 void SurfaceZPlane_to_hdf5_inner(hid_t group_id) const;   
 void SurfacePlane_to_hdf5_inner(hid_t group_id) const;    
 void SurfaceXCylinder_to_hdf5_inner(hid_t group_id) const;
 void SurfaceYCylinder_to_hdf5_inner(hid_t group_id) const;
 void SurfaceZCylinder_to_hdf5_inner(hid_t group_id) const;
 void SurfaceSphere_to_hdf5_inner(hid_t group_id) const;   
 void SurfaceXCone_to_hdf5_inner(hid_t group_id) const;    
 void SurfaceYCone_to_hdf5_inner(hid_t group_id) const;    
 void SurfaceZCone_to_hdf5_inner(hid_t group_id) const;    
 void SurfaceQuadric_to_hdf5_inner(hid_t group_id) const;  
};

//==============================================================================
//! A `Surface` representing a DAGMC-based surface in DAGMC.
//==============================================================================
/*
#ifdef DAGMC
class DAGSurface : public Surface
{
public:
  DAGSurface();

  double evaluate(Position r) const;
  double distance(Position r, Direction u, bool coincident) const;
  Direction normal(Position r) const;
  Direction reflect(Position r, Direction u, Particle* p) const;

  void to_hdf5(hid_t group_id) const;

  moab::DagMC* dagmc_ptr_; //!< Pointer to DagMC instance
  int32_t dag_index_;      //!< DagMC index of surface
};
#endif
*/

//==============================================================================
//! A `Surface` that supports periodic boundary conditions.
//!
//! Translational periodicity is supported for the `XPlane`, `YPlane`, `ZPlane`,
//! and `Plane` types.  Rotational periodicity is supported for
//! `XPlane`-`YPlane` pairs.
//==============================================================================

/*
class PeriodicSurface : public Surface
{
public:
  //int i_periodic_{C_NONE};    //!< Index of corresponding periodic surface

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
  BoundingBox bounding_box(bool pos_side) const;

  //double x0_;
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
  BoundingBox bounding_box(bool pos_side) const;

  //double y0_;
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
  BoundingBox bounding_box(bool pos_side) const;

  //double z0_;
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

  //double A_, B_, C_, D_;
};

//==============================================================================
//! A cylinder aligned along the x-axis.
//
//! The cylinder is described by the equation
//! \f$(y - y_0)^2 + (z - z_0)^2 - R^2 = 0\f$
//==============================================================================

class SurfaceXCylinder : public Surface
{
public:
  explicit SurfaceXCylinder(pugi::xml_node surf_node);
  double evaluate(Position r) const;
  double distance(Position r, Direction u, bool coincident) const;
  Direction normal(Position r) const;
  void to_hdf5_inner(hid_t group_id) const;
  BoundingBox bounding_box(bool pos_side) const;

  //double y0_, z0_, radius_;
};

//==============================================================================
//! A cylinder aligned along the y-axis.
//
//! The cylinder is described by the equation
//! \f$(x - x_0)^2 + (z - z_0)^2 - R^2 = 0\f$
//==============================================================================

class SurfaceYCylinder : public Surface
{
public:
  explicit SurfaceYCylinder(pugi::xml_node surf_node);
  double evaluate(Position r) const;
  double distance(Position r, Direction u, bool coincident) const;
  Direction normal(Position r) const;
  void to_hdf5_inner(hid_t group_id) const;
  BoundingBox bounding_box(bool pos_side) const;

  //double x0_, z0_, radius_;
};

//==============================================================================
//! A cylinder aligned along the z-axis.
//
//! The cylinder is described by the equation
//! \f$(x - x_0)^2 + (y - y_0)^2 - R^2 = 0\f$
//==============================================================================

class SurfaceZCylinder : public Surface
{
public:
  explicit SurfaceZCylinder(pugi::xml_node surf_node);
  double evaluate(Position r) const;
  double distance(Position r, Direction u, bool coincident) const;
  Direction normal(Position r) const;
  void to_hdf5_inner(hid_t group_id) const;
  BoundingBox bounding_box(bool pos_side) const;

  //double x0_, y0_, radius_;
};

//==============================================================================
//! A sphere.
//
//! The cylinder is described by the equation
//! \f$(x - x_0)^2 + (y - y_0)^2 + (z - z_0)^2 - R^2 = 0\f$
//==============================================================================

class SurfaceSphere : public Surface
{
public:
  explicit SurfaceSphere(pugi::xml_node surf_node);
  double evaluate(Position r) const;
  double distance(Position r, Direction u, bool coincident) const;
  Direction normal(Position r) const;
  void to_hdf5_inner(hid_t group_id) const;
  BoundingBox bounding_box(bool pos_side) const;

  //double x0_, y0_, z0_, radius_;
};

//==============================================================================
//! A cone aligned along the x-axis.
//
//! The cylinder is described by the equation
//! \f$(y - y_0)^2 + (z - z_0)^2 - R^2 (x - x_0)^2 = 0\f$
//==============================================================================

class SurfaceXCone : public Surface
{
public:
  explicit SurfaceXCone(pugi::xml_node surf_node);
  double evaluate(Position r) const;
  double distance(Position r, Direction u, bool coincident) const;
  Direction normal(Position r) const;
  void to_hdf5_inner(hid_t group_id) const;

  //double x0_, y0_, z0_, radius_sq_;
};

//==============================================================================
//! A cone aligned along the y-axis.
//
//! The cylinder is described by the equation
//! \f$(x - x_0)^2 + (z - z_0)^2 - R^2 (y - y_0)^2 = 0\f$
//==============================================================================

class SurfaceYCone : public Surface
{
public:
  explicit SurfaceYCone(pugi::xml_node surf_node);
  double evaluate(Position r) const;
  double distance(Position r, Direction u, bool coincident) const;
  Direction normal(Position r) const;
  void to_hdf5_inner(hid_t group_id) const;

  //double x0_, y0_, z0_, radius_sq_;
};

//==============================================================================
//! A cone aligned along the z-axis.
//
//! The cylinder is described by the equation
//! \f$(x - x_0)^2 + (y - y_0)^2 - R^2 (z - z_0)^2 = 0\f$
//==============================================================================

class SurfaceZCone : public Surface
{
public:
  explicit SurfaceZCone(pugi::xml_node surf_node);
  double evaluate(Position r) const;
  double distance(Position r, Direction u, bool coincident) const;
  Direction normal(Position r) const;
  void to_hdf5_inner(hid_t group_id) const;

  //double x0_, y0_, z0_, radius_sq_;
};

//==============================================================================
//! A general surface described by a quadratic equation.
//
//! \f$A x^2 + B y^2 + C z^2 + D x y + E y z + F x z + G x + H y + J z + K = 0\f$
//==============================================================================

class SurfaceQuadric : public Surface
{
public:
  explicit SurfaceQuadric(pugi::xml_node surf_node);
  double evaluate(Position r) const;
  double distance(Position r, Direction u, bool coincident) const;
  Direction normal(Position r) const;
  void to_hdf5_inner(hid_t group_id) const;

  // Ax^2 + By^2 + Cz^2 + Dxy + Eyz + Fxz + Gx + Hy + Jz + K = 0
  //double A_, B_, C_, D_, E_, F_, G_, H_, J_, K_;
};
*/

//==============================================================================
// Non-member functions
//==============================================================================

void read_surfaces(pugi::xml_node node);

void free_memory_surfaces();

} // namespace openmc
#endif // OPENMC_SURFACE_H
