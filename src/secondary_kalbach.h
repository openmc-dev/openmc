#ifndef OPENMC_SECONDARY_KALBACH_H
#define OPENMC_SECONDARY_KALBACH_H

#include <vector>

#include "hdf5.h"
#include "xtensor/xtensor.hpp"
#include "angle_energy.h"
#include "constants.h"
#include "endf.h"

namespace openmc {

class KalbachMann : public AngleEnergy {
public:
  explicit KalbachMann(hid_t group);
  void sample(double E_in, double& E_out, double& mu) const;
private:
  struct KMTable {
    int n_discrete;
    Interpolation interpolation;
    xt::xtensor<double, 1> e_out;
    xt::xtensor<double, 1> p;
    xt::xtensor<double, 1> c;
    xt::xtensor<double, 1> r;
    xt::xtensor<double, 1> a;
  };

  int n_region_;
  std::vector<int> breakpoints_;
  std::vector<Interpolation> interpolation_;
  std::vector<double> energy_;
  std::vector<KMTable> distribution_;
};

}

#endif // OPENMC_SECONDARY_KALBACH_H
