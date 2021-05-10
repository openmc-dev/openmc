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

  // Types of BoundaryCondition
  enum class BCType {
    Transmission,
    Vacuum,
    Reflective,
    White,
    TranslationalPeriodic,
    RotationalPeriodic
  };
  
  BoundaryCondition() {type_ = BCType::Transmission;}
  BoundaryCondition(BCType type);
  BoundaryCondition(BCType type, int i_surf, int j_surf);

  void init_TranslationalPeriodicBC();
  void init_RotationalPeriodicBC();

  //! Perform tracking operations for a particle that strikes the boundary.
  //! \param p The particle that struck the boundary.  This class is not meant
  //!   to directly modify anything about the particle, but it will do so
  //!   indirectly by calling the particle's appropriate cross_*_bc function.
  //! \param surf The specific surface on the boundary the particle struck.
  #pragma omp declare target
  void handle_particle(Particle& p, const Surface& surf) const;
  #pragma omp end declare target
void VacuumBC_handle_particle(Particle& p, const Surface& surf) const;
  #pragma omp declare target
void ReflectiveBC_handle_particle(Particle& p, const Surface& surf) const;
  #pragma omp end declare target
void WhiteBC_handle_particle(Particle& p, const Surface& surf) const;
void TranslationalPeriodicBC_handle_particle(Particle& p, const Surface& surf) const;
void RotationalPeriodicBC_handle_particle(Particle& p, const Surface& surf) const;

  //! Return a string classification of this BC.
  std::string type() const;

  BCType type_;
  int i_surf_;
  int j_surf_;
  Position translation_;
  double angle_;
};

/*
//==============================================================================
//! A BC that kills particles, indicating they left the problem.
//==============================================================================

class VacuumBC : public BoundaryCondition {
public:
  void
  handle_particle(Particle& p, const Surface& surf) const override;
};

//==============================================================================
//! A BC that returns particles via specular reflection.
//==============================================================================

class ReflectiveBC : public BoundaryCondition {
public:
  void
  handle_particle(Particle& p, const Surface& surf) const override;

};

//==============================================================================
//! A BC that returns particles via diffuse reflection.
//==============================================================================

class WhiteBC : public BoundaryCondition {
public:
  void
  handle_particle(Particle& p, const Surface& surf) const override;

};

//==============================================================================
//! A BC that moves particles to another part of the problem.
//==============================================================================

class PeriodicBC : public BoundaryCondition {
public:
  PeriodicBC(int i_surf, int j_surf)
    : i_surf_(i_surf), j_surf_(j_surf)
  {};

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

  void
  handle_particle(Particle& p, const Surface& surf) const override;

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

  void
  handle_particle(Particle& p, const Surface& surf) const override;

protected:
  //! Angle about the axis by which particle coordinates will be rotated
  double angle_;
};
*/

} // namespace openmc
#endif // OPENMC_BOUNDARY_CONDITION_H
