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

  BoundaryCondition() = default;
  BoundaryCondition(BCType type) : type_(type) {};
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
  void VacuumBC_handle_particle(Particle& p, const Surface& surf) const;
  void ReflectiveBC_handle_particle(Particle& p, const Surface& surf) const;
  void WhiteBC_handle_particle(Particle& p, const Surface& surf) const;
  #pragma omp end declare target
  void TranslationalPeriodicBC_handle_particle(Particle& p, const Surface& surf) const;
  void RotationalPeriodicBC_handle_particle(Particle& p, const Surface& surf) const;

  //! Return a string classification of this BC.
  std::string type() const;

  BCType type_ {BCType::Transmission};
  int i_surf_;
  int j_surf_;
  Position translation_;
  double angle_;
};

} // namespace openmc
#endif // OPENMC_BOUNDARY_CONDITION_H
