#ifndef OPENMC_DISTRIBUTION_SPATIAL_H
#define OPENMC_DISTRIBUTION_SPATIAL_H

#include "pugixml.hpp"

#include "openmc/distribution.h"
#include "openmc/position.h"

namespace openmc {

//==============================================================================
//! Probability density function for points in Euclidean space
//==============================================================================

class SpatialDistribution {
public:
  virtual ~SpatialDistribution() = default;

  //! Sample a position from the distribution
  virtual Position sample(uint64_t* seed) const = 0;
};

//==============================================================================
//! Distribution of points specified by independent distributions in x,y,z
//==============================================================================

class CartesianIndependent : public SpatialDistribution {
public:
  explicit CartesianIndependent(pugi::xml_node node);

  //! Sample a position from the distribution
  //! \param seed Pseudorandom number seed pointer
  //! \return Sampled position
  Position sample(uint64_t* seed) const;
private:
  UPtrDist x_; //!< Distribution of x coordinates
  UPtrDist y_; //!< Distribution of y coordinates
  UPtrDist z_; //!< Distribution of z coordinates
};

//==============================================================================
//! Distribution of points specified by spherical coordinates r,theta,phi
//==============================================================================

class SphericalIndependent : public SpatialDistribution {
public:
  explicit SphericalIndependent(pugi::xml_node node);

  //! Sample a position from the distribution
  //! \param seed Pseudorandom number seed pointer
  //! \return Sampled position
  Position sample(uint64_t* seed) const;
private:
  UPtrDist r_; //!< Distribution of r coordinates
  UPtrDist theta_; //!< Distribution of theta coordinates
  UPtrDist phi_; //!< Distribution of phi coordinates
  Position origin_; //!< Cartesian coordinates of the sphere center
};

//==============================================================================
//! Uniform distribution of points over a box
//==============================================================================

class SpatialBox : public SpatialDistribution {
public:
  explicit SpatialBox(pugi::xml_node node, bool fission=false);

  //! Sample a position from the distribution
  //! \param seed Pseudorandom number seed pointer
  //! \return Sampled position
  Position sample(uint64_t* seed) const;

  // Properties
  bool only_fissionable() const { return only_fissionable_; }
private:
  Position lower_left_; //!< Lower-left coordinates of box
  Position upper_right_; //!< Upper-right coordinates of box
  bool only_fissionable_ {false}; //!< Only accept sites in fissionable region?
};

//==============================================================================
//! Distribution at a single point
//==============================================================================

class SpatialPoint : public SpatialDistribution {
public:
  SpatialPoint() : r_{} { };
  SpatialPoint(Position r) : r_{r} { };
  explicit SpatialPoint(pugi::xml_node node);

  //! Sample a position from the distribution
  //! \param seed Pseudorandom number seed pointer
  //! \return Sampled position
  Position sample(uint64_t* seed) const;
private:
  Position r_; //!< Single position at which sites are generated
};

using UPtrSpace = std::unique_ptr<SpatialDistribution>;

} // namespace openmc

#endif // OPENMC_DISTRIBUTION_SPATIAL_H
