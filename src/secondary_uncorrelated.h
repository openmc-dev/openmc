#ifndef OPENMC_SECONDARY_UNCORRELATED_H
#define OPENMC_SECONDARY_UNCORRELATED_H

#include <memory>
#include <vector>

#include "hdf5.h"
#include "angle_energy.h"
#include "distribution_angle.h"
#include "distribution_energy.h"

namespace openmc {

class UncorrelatedAngleEnergy : public AngleEnergy {
public:
  explicit UncorrelatedAngleEnergy(hid_t group);
  void sample(double E_in, double& E_out, double& mu) const;
  double sampleMu(double E) const { return angle_.sample(E); }

  bool fission_ {false};
private:
  AngleDistribution angle_;
  std::unique_ptr<EnergyDistribution> energy_;
};

}

#endif // OPENMC_SECONDARY_UNCORRELATED_H
