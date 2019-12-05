//! \file distribution_angle.h
//! Angle distribution dependent on incident particle energy

#ifndef OPENMC_DISTRIBUTION_ANGLE_H
#define OPENMC_DISTRIBUTION_ANGLE_H

#include <vector> // for vector

#include "hdf5.h"

#include "openmc/distribution.h"

namespace openmc {

//==============================================================================
//! Angle distribution that depends on incident particle energy
//==============================================================================

class AngleDistribution {
public:
  AngleDistribution() = default;
  explicit AngleDistribution(hid_t group);

  //! Sample an angle given an incident particle energy
  //! \param[in] E Particle energy in [eV]
  //! \param[inout] seed pseudorandom number seed pointer
  //! \return Cosine of the angle in the range [-1,1]
  double sample(double E, uint64_t* seed) const;

  //! Determine whether angle distribution is empty
  //! \return Whether distribution is empty
  bool empty() const { return energy_.empty(); }

private:
  std::vector<double> energy_;
  std::vector<UPtrDist> distribution_;
};

} // namespace openmc

#endif // OPENMC_DISTRIBUTION_ANGLE_H
