#ifndef OPENMC_SECONDARY_CORRELATED_H
#define OPENMC_SECONDARY_CORRELATED_H

#include <vector>

#include "hdf5.h"
#include "xtensor/xtensor.hpp"
#include "angle_energy.h"
#include "endf.h"
#include "distribution.h"

namespace openmc {

class CorrelatedAngleEnergy : public AngleEnergy {
public:
  explicit CorrelatedAngleEnergy(hid_t group);
  void sample(double E_in, double& E_out, double& mu) const;
private:
  struct CorrTable {
    int n_discrete;
    Interpolation interpolation;
    xt::xtensor<double, 1> e_out;
    xt::xtensor<double, 1> p;
    xt::xtensor<double, 1> c;
    std::vector<UPtrDist> angle;
  };

  int n_region_;
  std::vector<int> breakpoints_;
  std::vector<Interpolation> interpolation_;
  std::vector<double> energy_;
  std::vector<CorrTable> distribution_;
};

}

#endif // OPENMC_SECONDARY_CORRELATED_H
