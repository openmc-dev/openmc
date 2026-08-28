#ifndef OPENMC_BOUNDARY_CONDITION_H
#define OPENMC_BOUNDARY_CONDITION_H

#include "openmc/hdf5_interface.h"
#include "openmc/particle.h"
#include "openmc/position.h"
#include <fmt/core.h>

namespace openmc {

// Forward declare some types used in function arguments.
class Particle;
class RandomRay;
class Surface;

//==============================================================================
//! A class that tells particles what to do after they strike an outer boundary.
//==============================================================================

class BoundaryCondition {
public:
  virtual ~BoundaryCondition() = default;

  //! Perform tracking operations for a particle that strikes the boundary.
  //! \param p The particle that struck the boundary.  This class is not meant
  //!   to directly modify anything about the particle, but it will do so
  //!   indirectly by calling the particle's appropriate cross_*_bc function.
  //! \param surf The specific surface on the boundary the particle struck.
  virtual void handle_particle(Particle& p, const Surface& surf) const = 0;

  //! Modify the incident particle's weight according to the boundary's albedo.
  //! \param p The particle that struck the boundary.  This function calculates
  //!   the reduction in the incident particle's weight as it interacts
  //!   with a boundary. The lost weight is tallied before the remaining weight
  //!   is reassigned to the incident particle. Implementations of the
  //!   handle_particle function typically call this method in its body.
  //! \param surf The specific surface on the boundary the particle struck.
  void handle_albedo(Particle& p, const Surface& surf) const
  {
    if (!has_albedo())
      return;
    double initial_wgt = p.wgt();
    // Treat the lost weight fraction as leakage, similar to VacuumBC.
    // This ensures the lost weight is tallied properly.
    p.wgt() *= (1.0 - albedo_);
    p.cross_vacuum_bc(surf);
    p.wgt() = initial_wgt * albedo_;
  };

  //! Return a string classification of this BC.
  virtual std::string type() const = 0;

  //! Write albedo data of this BC to hdf5.
  void to_hdf5(hid_t surf_group) const
  {
    if (has_albedo()) {
      write_string(surf_group, "albedo", fmt::format("{}", albedo_), false);
    }
  };

  //! Set albedo of this BC.
  void set_albedo(double albedo) { albedo_ = albedo; }

  //! Return if this BC has an albedo.
  bool has_albedo() const { return (albedo_ > 0.0); }

  //! Return albedo value.
  double albedo() const { return albedo_; }

private:
  double albedo_ = -1.0;
};

//==============================================================================
//! A BC that kills particles, indicating they left the problem.
//==============================================================================

class VacuumBC : public BoundaryCondition {
public:
  void handle_particle(Particle& p, const Surface& surf) const override;

  std::string type() const override { return "vacuum"; }
};

//==============================================================================
//! A BC that transforms particle direction and/or location.
//==============================================================================

class TransformingBC : public BoundaryCondition {
public:
  //! Apply the geometric transformation.
  //! \param[in] p Particle
  //! \param[in] surf Current surface
  //! \param[out] new_r Transformed position
  //! \param[out] new_u Transformed direction
  //! \param[out] new_surface New surface index
  virtual void transform(Particle& p, const Surface& surf, Position& new_r,
    Direction& new_u, int& new_surface) const = 0;

  //! Handle the crossing of the boundary condition by the particle
  //! \param[inout] p Particle
  //! \param[in] surf Current surface
  //! \param[in] new_r Transformed position
  //! \param[in] new_u Transformed direction
  //! \param[in] new_surface New surface index
  virtual void cross_bc(Particle& p, const Surface& surf, const Position& new_r,
    const Direction& new_u, int new_surface) const = 0;

  void handle_particle(Particle& p, const Surface& surf) const override
  {
    Position new_r;
    Direction new_u;
    int new_surface;
    transform(p, surf, new_r, new_u, new_surface);
    BoundaryCondition::handle_albedo(p, surf);
    cross_bc(p, surf, new_r, new_u, new_surface);
  }
};

//==============================================================================
//! A BC that returns particles via specular reflection.
//==============================================================================

class ReflectiveBC : public TransformingBC {
public:
  void transform(Particle& p, const Surface& surf, Position& new_r,
    Direction& new_u, int& new_surface) const override;

  void cross_bc(Particle& p, const Surface& surf, const Position& new_r,
    const Direction& new_u, int new_surface) const override;

  std::string type() const override { return "reflective"; }
};

//==============================================================================
//! A BC that returns particles via diffuse reflection.
//==============================================================================

class WhiteBC : public TransformingBC {
public:
  void transform(Particle& p, const Surface& surf, Position& new_r,
    Direction& new_u, int& new_surface) const override;

  void cross_bc(Particle& p, const Surface& surf, const Position& new_r,
    const Direction& new_u, int new_surface) const override;

  std::string type() const override { return "white"; }
};

//==============================================================================
//! A BC that moves particles to another part of the problem.
//==============================================================================

class PeriodicBC : public TransformingBC {
public:
  PeriodicBC(int i_surf, int j_surf) : i_surf_(i_surf), j_surf_(j_surf) {};

  std::string type() const override { return "periodic"; }

  int i_surf() const { return i_surf_; }

  int j_surf() const { return j_surf_; }

protected:
  int i_surf_;
  int j_surf_;
};

//==============================================================================
//! A BC that moves particles to another part of the problem without rotation.
//==============================================================================

class TranslationalPeriodicBC : public PeriodicBC {
public:
  TranslationalPeriodicBC(int i_surf, int j_surf);

  void transform(Particle& p, const Surface& surf, Position& new_r,
    Direction& new_u, int& new_surface) const override;

  void cross_bc(Particle& p, const Surface& surf, const Position& new_r,
    const Direction& new_u, int new_surface) const override;

protected:
  //! Vector along which incident particles will be moved
  Position translation_;
};

//==============================================================================
//! A BC that rotates particles about a global axis.
//
//! Only rotations about the x, y, and z axes are supported.
//==============================================================================

class RotationalPeriodicBC : public PeriodicBC {
public:
  enum PeriodicAxis { x, y, z };
  RotationalPeriodicBC(int i_surf, int j_surf, PeriodicAxis axis);

  void transform(Particle& p, const Surface& surf, Position& new_r,
    Direction& new_u, int& new_surface) const override;

  void cross_bc(Particle& p, const Surface& surf, const Position& new_r,
    const Direction& new_u, int new_surface) const override;

protected:
  //! Angle about the axis by which particle coordinates will be rotated
  double angle_;
  //! Do we need to flip surfaces senses when applying the transformation?
  bool flip_sense_;
  //! Ensure that choice of axes is right handed. axis_1_idx_ corresponds to the
  //! independent axis and axis_2_idx_ corresponds to the dependent axis in the
  //! 2D plane  perpendicular to the planes' axis of rotation
  int zero_axis_idx_;
  int axis_1_idx_;
  int axis_2_idx_;
};

} // namespace openmc
#endif // OPENMC_BOUNDARY_CONDITION_H
