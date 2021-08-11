#ifndef OPENMC_BOUNDARY_CONDITION_H
#define OPENMC_BOUNDARY_CONDITION_H

#include "openmc/position.h"

namespace openmc {

// Forward declare some types used in function arguments.
class Particle;
class Surface;

//==============================================================================
//! A class that tells particles what to do after they strike an outer boundary.
//==============================================================================

class BoundaryCondition {
public:
  //! Perform tracking operations for a particle that strikes the boundary.
  //! \param p The particle that struck the boundary.  This class is not meant
  //!   to directly modify anything about the particle, but it will do so
  //!   indirectly by calling the particle's appropriate cross_*_bc function.
  //! \param surf The specific surface on the boundary the particle struck.
  virtual void handle_particle(Particle& p, const Surface& surf) const = 0;

  //! Return a string classification of this BC.
  virtual std::string type() const = 0;
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
//! A BC that returns particles via specular reflection.
//==============================================================================

class ReflectiveBC : public BoundaryCondition {
public:
  void handle_particle(Particle& p, const Surface& surf) const override;

  std::string type() const override { return "reflective"; }
};

//==============================================================================
//! A BC that returns particles via diffuse reflection.
//==============================================================================

class WhiteBC : public BoundaryCondition {
public:
  void handle_particle(Particle& p, const Surface& surf) const override;

  std::string type() const override { return "white"; }
};

//==============================================================================
//! A BC that moves particles to another part of the problem.
//==============================================================================

class PeriodicBC : public BoundaryCondition {
public:
  PeriodicBC(int i_surf, int j_surf) : i_surf_(i_surf), j_surf_(j_surf) {};

  std::string type() const override { return "periodic"; }

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

  void handle_particle(Particle& p, const Surface& surf) const override;

protected:
  //! Vector along which incident particles will be moved
  Position translation_;
};

//==============================================================================
//! A BC that rotates particles about a global axis.
//
//! Currently only rotations about the z-axis are supported.
//==============================================================================

class RotationalPeriodicBC : public PeriodicBC {
public:
  RotationalPeriodicBC(int i_surf, int j_surf);

  void handle_particle(Particle& p, const Surface& surf) const override;

protected:
  //! Angle about the axis by which particle coordinates will be rotated
  double angle_;
};

} // namespace openmc
#endif // OPENMC_BOUNDARY_CONDITION_H
