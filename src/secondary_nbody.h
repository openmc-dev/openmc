#ifndef OPENMC_SECONDARY_NBODY_H
#define OPENMC_SECONDARY_NBODY_H

#include "hdf5.h"

#include "angle_energy.h"

namespace openmc {

class NBodyPhaseSpace : public AngleEnergy {
public:
  explicit NBodyPhaseSpace(hid_t group);
  void sample(double E_in, double& E_out, double& mu) const;
private:
  int n_bodies_;
  double mass_ratio_;
  double A_;
  double Q_;
};

}

#endif // OPENMC_SECONDARY_NBODY_H
