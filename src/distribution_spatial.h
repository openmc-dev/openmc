#ifndef OPENMC_DISTRIBTUION_SPATIAL_H
#define OPENMC_DISTRIBUTION_SPATIAL_H

#include "pugixml.hpp"

#include "distribution.h"
#include "position.h"

namespace openmc {

//==============================================================================
//! Probability density function for points in Euclidean space
//==============================================================================

class SpatialDistribution {
public:
  virtual Position sample() const = 0;
  virtual ~SpatialDistribution() = default;
};

//==============================================================================
//! Distribution of points specified by independent distributions in x,y,z
//==============================================================================

class CartesianIndependent : public SpatialDistribution {
public:
  explicit CartesianIndependent(pugi::xml_node node);

  //! Sample a position from the distribution
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
  Position sample() const;
private:
  Position lower_left_;
  Position upper_right_;
  bool only_fissionable {false};
};

//==============================================================================
//! Distribution at a single point
//==============================================================================

class SpatialPoint : public SpatialDistribution {
public:
  explicit SpatialPoint(pugi::xml_node node);

  //! Sample a position from the distribution
  Position sample() const;
private:
  Position r_;
};

} // namespace openmc

#endif // OPENMC_DISTRIBUTION_SPATIAL_H
