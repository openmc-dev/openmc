#ifndef DISTRIBUTION_MULTI_H
#define DISTRIBUTION_MULTI_H

#include "openmc/memory.h"

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
  UnitSphereDistribution() {};
  explicit UnitSphereDistribution(Direction u) : u_ref_ {u} {};
  explicit UnitSphereDistribution(pugi::xml_node node);
  virtual ~UnitSphereDistribution() = default;

  static unique_ptr<UnitSphereDistribution> create(pugi::xml_node node);

  //! Sample a direction from the distribution
  //! \param seed Pseudorandom number seed pointer
  //! \return (sampled Direction, sample weight)
  virtual std::pair<Direction, double> sample(uint64_t* seed) const = 0;

  Direction u_ref_ {0.0, 0.0, 1.0}; //!< reference direction
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
  //! \return (sampled Direction, sample weight)
  std::pair<Direction, double> sample(uint64_t* seed) const override;

  //! Sample a direction and return evaluation of the PDF for biased sampling.
  //! Note that bias distributions are intended to return unit-weight samples.
  //! \param seed Pseudorandom number seed points
  //! \return (sampled Direction, value of the PDF at this Direction)
  std::pair<Direction, double> sample_as_bias(uint64_t* seed) const;

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

Direction isotropic_direction(uint64_t* seed);

class Isotropic : public UnitSphereDistribution {
public:
  Isotropic() {};
  explicit Isotropic(pugi::xml_node node);

  //! Sample a direction from the distribution
  //! \param seed Pseudorandom number seed pointer
  //! \return (sampled Direction, sample weight)
  std::pair<Direction, double> sample(uint64_t* seed) const override;

  // Set or get bias distribution
  void set_bias(std::unique_ptr<PolarAzimuthal> bias)
  {
    bias_ = std::move(bias);
  }

  const PolarAzimuthal* bias() const { return bias_.get(); }

protected:
  // Biasing distribution
  unique_ptr<PolarAzimuthal> bias_;
};

//==============================================================================
//! Monodirectional distribution
//==============================================================================

class Monodirectional : public UnitSphereDistribution {
public:
  Monodirectional(Direction u) : UnitSphereDistribution {u} {};
  explicit Monodirectional(pugi::xml_node node)
    : UnitSphereDistribution {node} {};

  //! Sample a direction from the distribution
  //! \param seed Pseudorandom number seed pointer
  //! \return (sampled Direction, sample weight)
  std::pair<Direction, double> sample(uint64_t* seed) const override;
};

using UPtrAngle = unique_ptr<UnitSphereDistribution>;

} // namespace openmc

#endif // DISTRIBUTION_MULTI_H
