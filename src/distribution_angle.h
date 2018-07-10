#ifndef OPENMC_DISTRIBUTION_ANGLE_H
#define OPENMC_DISTRIBUTION_ANGLE_H

#include <vector> // for vector

#include "distribution.h"
#include "hdf5.h"

namespace openmc {

class AngleDistribution {
public:
  AngleDistribution() = default;
  explicit AngleDistribution(hid_t group);
  double sample(double E) const;

  bool empty() const { return energy_.empty(); }
private:
  std::vector<double> energy_;
  std::vector<UPtrDist> distribution_;
};

}

#endif // OPENMC_DISTRIBUTION_ANGLE_H
