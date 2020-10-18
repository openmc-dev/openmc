#ifndef OPENMC_BOUNDARY_CONDITION_H
#define OPENMC_BOUNDARY_CONDITION_H

#include "openmc/position.h"

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

//class BoundaryCondition;
//
//namespace model {
//  extern std::vector<std::unique_ptr<BoundaryCondition>> boundary_conditions;
//} // namespace model

class Particle;
class Surface;

class BoundaryCondition
{
public:
  virtual void
  handle_particle(Particle& p, const Surface& surf) const = 0;
};

class VacuumBC : public BoundaryCondition
{
public:
  void
  handle_particle(Particle& p, const Surface& surf) const override;
};

class ReflectiveBC : public BoundaryCondition
{
public:
  void
  handle_particle(Particle& p, const Surface& surf) const override;
};

class WhiteBC : public BoundaryCondition
{
public:
  void
  handle_particle(Particle& p, const Surface& surf) const override;
};

class PeriodicBC : public BoundaryCondition
{
public:
  PeriodicBC(int i_surf, int j_surf)
    : i_surf_(i_surf), j_surf_(j_surf)
  {};

protected:
  int i_surf_;
  int j_surf_;
};

class TranslationalPeriodicBC : public PeriodicBC
{
public:
  TranslationalPeriodicBC(int i_surf, int j_surf);

  void
  handle_particle(Particle& p, const Surface& surf) const override;

protected:
  Position translation_;
};

class RotationalPeriodicBC : public PeriodicBC
{
public:
  RotationalPeriodicBC(int i_surf, int j_surf);

  void
  handle_particle(Particle& p, const Surface& surf) const override;
};

} // namespace openmc
#endif // OPENMC_BOUNDARY_CONDITION_H
