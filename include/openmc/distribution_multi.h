#ifndef DISTRIBUTION_MULTI_H
#define DISTRIBUTION_MULTI_H

#include <memory>

#include "openmc/distribution.h"
#include "openmc/position.h"

namespace openmc {

//==============================================================================
//! Probability density function for points on the unit sphere. Extensions of
//! this type are used to sample angular distributions for starting sources
//==============================================================================

class UnitSphereDistribution {
public:
  UnitSphereDistribution() { };
  explicit UnitSphereDistribution(Direction u) : u_ref{u} { };
  virtual ~UnitSphereDistribution() = default;

  //! Sample a direction from the distribution
  //! \return Direction sampled
  virtual Direction sample() const = 0;

  Direction u_ref {0.0, 0.0, 1.0};  //!< reference direction
};

//==============================================================================
//! Explicit distribution of polar and azimuthal angles
//==============================================================================

class PolarAzimuthal : public UnitSphereDistribution {
public:
  PolarAzimuthal(Direction u, UPtrDist mu, UPtrDist phi);

  //! Sample a direction from the distribution
  //! \return Direction sampled
  Direction sample() const;
private:
  UPtrDist mu_;  //!< Distribution of polar angle
  UPtrDist phi_; //!< Distribution of azimuthal angle
};

//==============================================================================
//! Uniform distribution on the unit sphere
//==============================================================================

class Isotropic : public UnitSphereDistribution {
public:
  Isotropic() { };

  //! Sample a direction from the distribution
  //! \return Sampled direction
  Direction sample() const;
};

//==============================================================================
//! Monodirectional distribution
//==============================================================================

class Monodirectional : public UnitSphereDistribution {
public:
  Monodirectional(Direction u) : UnitSphereDistribution{u} { };

  //! Sample a direction from the distribution
  //! \return Sampled direction
  Direction sample() const;
};

} // namespace openmc

#endif // DISTRIBUTION_MULTI_H
