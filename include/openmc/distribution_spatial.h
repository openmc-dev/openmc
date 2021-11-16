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
  HD virtual Position sample(uint64_t* seed) const = 0;
  HD virtual bool only_fissionable() const { return false; }
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
  HD virtual Position sample(uint64_t* seed) const override;

  // Observer pointers
  Distribution* x() const { return x_.get(); }
  Distribution* y() const { return y_.get(); }
  Distribution* z() const { return z_.get(); }
private:
  unique_ptr<Distribution> x_; //!< Distribution of x coordinates
  unique_ptr<Distribution> y_; //!< Distribution of y coordinates
  unique_ptr<Distribution> z_; //!< Distribution of z coordinates
};

//==============================================================================
//! Distribution of points specified by cylindrical coordinates r,phi,z
//==============================================================================

class CylindricalIndependent : public SpatialDistribution {
public:
  explicit CylindricalIndependent(pugi::xml_node node);

  //! Sample a position from the distribution
  //! \param seed Pseudorandom number seed pointer
  //! \return Sampled position
  HD virtual Position sample(uint64_t* seed) const;

  Distribution* r() const { return r_.get(); }
  Distribution* phi() const { return phi_.get(); }
  Distribution* z() const { return z_.get(); }
  Position origin() const { return origin_; }
private:
  unique_ptr<Distribution> r_;   //!< Distribution of r coordinates
  unique_ptr<Distribution> phi_; //!< Distribution of phi coordinates
  unique_ptr<Distribution> z_;   //!< Distribution of z coordinates
  Position origin_; //!< Cartesian coordinates of the cylinder center
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
  HD virtual Position sample(uint64_t* seed) const;

  Distribution* r() const { return r_.get(); }
  Distribution* theta() const { return theta_.get(); }
  Distribution* phi() const { return phi_.get(); }
  Position origin () const { return origin_; }
private:
  unique_ptr<Distribution> r_;     //!< Distribution of r coordinates
  unique_ptr<Distribution> theta_; //!< Distribution of theta coordinates
  unique_ptr<Distribution> phi_;   //!< Distribution of phi coordinates
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
  HD virtual Position sample(uint64_t* seed) const override;

  // Properties
  HD virtual bool only_fissionable() const override
  {
    return only_fissionable_;
  }
  Position lower_left() const { return lower_left_; }
  Position upper_right() const { return upper_right_; }
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
  HD virtual Position sample(uint64_t* seed) const override;

  Position r() const { return r_; }
private:
  Position r_; //!< Single position at which sites are generated
};

using UPtrSpace = unique_ptr<SpatialDistribution>;

} // namespace openmc

#endif // OPENMC_DISTRIBUTION_SPATIAL_H
