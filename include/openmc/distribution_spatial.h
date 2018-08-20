#ifndef OPENMC_DISTRIBTUION_SPATIAL_H
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
  virtual Position sample() const = 0;
};

//==============================================================================
//! Distribution of points specified by independent distributions in x,y,z
//==============================================================================

class CartesianIndependent : public SpatialDistribution {
public:
  explicit CartesianIndependent(pugi::xml_node node);

  //! Sample a position from the distribution
  //! \return Sampled position
  Position sample() const;
private:
  UPtrDist x_; //!< Distribution of x coordinates
  UPtrDist y_; //!< Distribution of y coordinates
  UPtrDist z_; //!< Distribution of z coordinates
};

//==============================================================================
//! Uniform distribution of points over a box
//==============================================================================

class SpatialBox : public SpatialDistribution {
public:
  explicit SpatialBox(pugi::xml_node node);

  //! Sample a position from the distribution
  //! \return Sampled position
  Position sample() const;
private:
  Position lower_left_; //!< Lower-left coordinates of box
  Position upper_right_; //!< Upper-right coordinates of box
  bool only_fissionable {false}; //!< Only accept sites in fissionable region?
};

//==============================================================================
//! Distribution at a single point
//==============================================================================

class SpatialPoint : public SpatialDistribution {
public:
  explicit SpatialPoint(pugi::xml_node node);

  //! Sample a position from the distribution
  //! \return Sampled position
  Position sample() const;
private:
  Position r_; //!< Single position at which sites are generated
};

} // namespace openmc

#endif // OPENMC_DISTRIBUTION_SPATIAL_H
