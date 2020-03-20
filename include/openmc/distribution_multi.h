#ifndef DISTRIBUTION_MULTI_H
#define DISTRIBUTION_MULTI_H

#include <memory>

#include "pugixml.hpp"

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
  explicit UnitSphereDistribution(Direction u) : u_ref_{u} { };
  explicit UnitSphereDistribution(pugi::xml_node node);
  virtual ~UnitSphereDistribution() = default;

  //! Sample a direction from the distribution
  //! \param seed Pseudorandom number seed pointer
  //! \return Direction sampled
  virtual Direction sample(uint64_t* seed) const = 0;

  Direction u_ref_ {0.0, 0.0, 1.0};  //!< reference direction
};

//==============================================================================
//! Explicit distribution of polar and azimuthal angles
//==============================================================================

class PolarAzimuthal : public UnitSphereDistribution {
public:
  PolarAzimuthal(Direction u, UPtrDist mu, UPtrDist phi);
  explicit PolarAzimuthal(pugi::xml_node node);

  //! Sample a direction from the distribution
  //! \param seed Pseudorandom number seed pointer
  //! \return Direction sampled
  Direction sample(uint64_t* seed) const;

  // Observing pointers
  Distribution* mu() const { return mu_.get(); }
  Distribution* phi() const { return phi_.get(); }

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
  //! \param seed Pseudorandom number seed pointer
  //! \return Sampled direction
  Direction sample(uint64_t* seed) const;
};

//==============================================================================
//! Monodirectional distribution
//==============================================================================

class Monodirectional : public UnitSphereDistribution {
public:
  Monodirectional(Direction u) : UnitSphereDistribution{u} { };
  explicit Monodirectional(pugi::xml_node node) : UnitSphereDistribution{node} { };

  //! Sample a direction from the distribution
  //! \param seed Pseudorandom number seed pointer
  //! \return Sampled direction
  Direction sample(uint64_t* seed) const;
};

using UPtrAngle = std::unique_ptr<UnitSphereDistribution>;

} // namespace openmc

#endif // DISTRIBUTION_MULTI_H
