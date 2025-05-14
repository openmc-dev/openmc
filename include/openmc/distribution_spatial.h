#ifndef OPENMC_DISTRIBUTION_SPATIAL_H
#define OPENMC_DISTRIBUTION_SPATIAL_H

#include "pugixml.hpp"

#include "openmc/distribution.h"
#include "openmc/mesh.h"
#include "openmc/position.h"
#include "openmc/span.h"

namespace openmc {

//==============================================================================
//! Probability density function for points in Euclidean space
//==============================================================================

class SpatialDistribution {
public:
  virtual ~SpatialDistribution() = default;

  //! Sample a position from the distribution
  virtual Position sample(uint64_t* seed) const = 0;

  static unique_ptr<SpatialDistribution> create(pugi::xml_node node);
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
  Position sample(uint64_t* seed) const override;

  // Observer pointers
  Distribution* x() const { return x_.get(); }
  Distribution* y() const { return y_.get(); }
  Distribution* z() const { return z_.get(); }

private:
  UPtrDist x_; //!< Distribution of x coordinates
  UPtrDist y_; //!< Distribution of y coordinates
  UPtrDist z_; //!< Distribution of z coordinates
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
  Position sample(uint64_t* seed) const override;

  Distribution* r() const { return r_.get(); }
  Distribution* phi() const { return phi_.get(); }
  Distribution* z() const { return z_.get(); }
  Position origin() const { return origin_; }

private:
  UPtrDist r_;      //!< Distribution of r coordinates
  UPtrDist phi_;    //!< Distribution of phi coordinates
  UPtrDist z_;      //!< Distribution of z coordinates
  Position origin_; //!< Cartesian coordinates of the cylinder center
};

//==============================================================================
//! Distribution of points specified by spherical coordinates r,cos_theta,phi
//==============================================================================

class SphericalIndependent : public SpatialDistribution {
public:
  explicit SphericalIndependent(pugi::xml_node node);

  //! Sample a position from the distribution
  //! \param seed Pseudorandom number seed pointer
  //! \return Sampled position
  Position sample(uint64_t* seed) const override;

  Distribution* r() const { return r_.get(); }
  Distribution* cos_theta() const { return cos_theta_.get(); }
  Distribution* phi() const { return phi_.get(); }
  Position origin() const { return origin_; }

private:
  UPtrDist r_;         //!< Distribution of r coordinates
  UPtrDist cos_theta_; //!< Distribution of cos_theta coordinates
  UPtrDist phi_;       //!< Distribution of phi coordinates
  Position origin_;    //!< Cartesian coordinates of the sphere center
};

//==============================================================================
//! Distribution of points within a mesh
//==============================================================================

class MeshSpatial : public SpatialDistribution {
public:
  explicit MeshSpatial(pugi::xml_node node);
  explicit MeshSpatial(int32_t mesh_id, span<const double> strengths);

  //! Sample a position from the distribution
  //! \param seed Pseudorandom number seed pointer
  //! \return Sampled position
  Position sample(uint64_t* seed) const override;

  //! Sample the mesh for an element and position within that element
  //! \param seed Pseudorandom number seed pointer
  //! \return Sampled element index and position within that element
  std::pair<int32_t, Position> sample_mesh(uint64_t* seed) const;

  //! Sample a mesh element
  //! \param seed Pseudorandom number seed pointer
  //! \return Sampled element index
  int32_t sample_element_index(uint64_t* seed) const;

  //! For unstructured meshes, ensure that elements are all linear tetrahedra
  void check_element_types() const;

  // Accessors
  const Mesh* mesh() const { return model::meshes.at(mesh_idx_).get(); }
  int32_t n_sources() const { return this->mesh()->n_bins(); }

  double total_strength() { return this->elem_idx_dist_.integral(); }

private:
  int32_t mesh_idx_ {C_NONE};
  DiscreteIndex elem_idx_dist_; //!< Distribution of
                                //!< mesh element indices
};

//==============================================================================
//! Distribution of points
//==============================================================================

class PointCloud : public SpatialDistribution {
public:
  explicit PointCloud(pugi::xml_node node);
  explicit PointCloud(
    std::vector<Position> point_cloud, span<const double> strengths);

  //! Sample a position from the distribution
  //! \param seed Pseudorandom number seed pointer
  //! \return Sampled position
  Position sample(uint64_t* seed) const override;

private:
  std::vector<Position> point_cloud_;
  DiscreteIndex point_idx_dist_; //!< Distribution of Position indices
};

//==============================================================================
//! Uniform distribution of points over a box
//==============================================================================

class SpatialBox : public SpatialDistribution {
public:
  explicit SpatialBox(pugi::xml_node node, bool fission = false);

  //! Sample a position from the distribution
  //! \param seed Pseudorandom number seed pointer
  //! \return Sampled position
  Position sample(uint64_t* seed) const override;

  // Properties
  bool only_fissionable() const { return only_fissionable_; }
  Position lower_left() const { return lower_left_; }
  Position upper_right() const { return upper_right_; }

private:
  Position lower_left_;           //!< Lower-left coordinates of box
  Position upper_right_;          //!< Upper-right coordinates of box
  bool only_fissionable_ {false}; //!< Only accept sites in fissionable region?
};

//==============================================================================
//! Uniform distribution of points over a ball
//==============================================================================

class SpatialBall : public SpatialDistribution {
public:
  explicit SpatialBall(pugi::xml_node node);

  //! Sample a position from the distribution
  //! \param seed Pseudorandom number seed pointer
  //! \return Sampled position
  Position sample(uint64_t* seed) const override;

  // Properties
  Position origin() const { return origin_; }
  double radius() const { return radius_; }

private:
  Position origin_; //!< Origin coordinates of ball
  double radius_;   //!< Radius of ball
};

//==============================================================================
//! Distribution at a single point
//==============================================================================

class SpatialPoint : public SpatialDistribution {
public:
  SpatialPoint() : r_ {} {};
  SpatialPoint(Position r) : r_ {r} {};
  explicit SpatialPoint(pugi::xml_node node);

  //! Sample a position from the distribution
  //! \param seed Pseudorandom number seed pointer
  //! \return Sampled position
  Position sample(uint64_t* seed) const override;

  Position r() const { return r_; }

private:
  Position r_; //!< Single position at which sites are generated
};

using UPtrSpace = unique_ptr<SpatialDistribution>;

} // namespace openmc

#endif // OPENMC_DISTRIBUTION_SPATIAL_H
