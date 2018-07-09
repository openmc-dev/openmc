#ifndef OPENMC_DISTRIBUTION_ANGLE_H
#define OPENMC_DISTRIBUTION_ANGLE_H

#include <memory> // for unique_ptr
#include <vector> // for vector

#include "distribution.h"
#include "hdf5.h"

namespace openmc {

class AngleDistribution {
public:
  explicit AngleDistribution(hid_t group);
  double sample(double E) const;
private:
  std::vector<double> energy_;
  std::vector<UPtrDist> distribution_;
};

}

#endif // OPENMC_DISTRIBUTION_ANGLE_H
