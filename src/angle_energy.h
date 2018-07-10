#ifndef OPENMC_ANGLE_ENERGY_H
#define OPENMC_ANGLE_ENERGY_H

namespace openmc {

class AngleEnergy {
public:
  virtual void sample(double E_in, double& E_out, double& mu) const = 0;
  virtual ~AngleEnergy() = default;
};

}

#endif // OPENMC_ANGLE_ENERGY_H
