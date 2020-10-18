#ifndef OPENMC_BOUNDARY_CONDITION_H
#define OPENMC_BOUNDARY_CONDITION_H

#include "openmc/position.h"

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

class Particle;
class Surface;

class BoundaryCondition
{
public:
  virtual void
  handle_particle(Particle& p, const Surface& surf) const = 0;

  virtual std::string type() const = 0;
};

class VacuumBC : public BoundaryCondition
{
public:
  void
  handle_particle(Particle& p, const Surface& surf) const override;

  std::string type() const override {return "vacuum";}
};

class ReflectiveBC : public BoundaryCondition
{
public:
  void
  handle_particle(Particle& p, const Surface& surf) const override;

  std::string type() const override {return "reflective";}
};

class WhiteBC : public BoundaryCondition
{
public:
  void
  handle_particle(Particle& p, const Surface& surf) const override;

  std::string type() const override {return "white";}
};

class PeriodicBC : public BoundaryCondition
{
public:
  PeriodicBC(int i_surf, int j_surf)
    : i_surf_(i_surf), j_surf_(j_surf)
  {};

  std::string type() const override {return "periodic";}

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

protected:
  double angle_;
};

} // namespace openmc
#endif // OPENMC_BOUNDARY_CONDITION_H
