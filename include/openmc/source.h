//! \file source.h
//! \brief External source distributions

#ifndef OPENMC_SOURCE_H
#define OPENMC_SOURCE_H

#include <memory>
#include <vector>

#include "pugixml.hpp"

#include "openmc/distribution_multi.h"
#include "openmc/distribution_spatial.h"
#include "openmc/particle.h"

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

class SourceDistribution;

namespace model {

extern std::vector<std::unique_ptr<SourceDistribution>> external_sources;

} // namespace model

//==============================================================================
//! External source distribution
//==============================================================================

class SourceDistribution {
public:
  virtual ~SourceDistribution() = default;

  // Methods that must be implemented
  virtual Particle::Bank sample(uint64_t* seed) = 0;

  // Methods that can be overridden
  virtual double strength() const { return 1.0; }
};

class IndependentSourceDistribution : public SourceDistribution {
public:
  // Constructors
  IndependentSourceDistribution(UPtrSpace space, UPtrAngle angle, UPtrDist energy);
  explicit IndependentSourceDistribution(pugi::xml_node node);

  //! Sample from the external source distribution
  //! \param[inout] seed Pseudorandom seed pointer
  //! \return Sampled site
  Particle::Bank sample(uint64_t* seed) override;

  // Properties
  Particle::Type particle_type() const { return particle_; }
  double strength() const override { return strength_; }

  // Make observing pointers available
  SpatialDistribution* space() const { return space_.get(); }
  UnitSphereDistribution* angle() const { return angle_.get(); }
  Distribution* energy() const { return energy_.get(); }

private:
  Particle::Type particle_ {Particle::Type::neutron}; //!< Type of particle emitted
  double strength_ {1.0}; //!< Source strength
  UPtrSpace space_; //!< Spatial distribution
  UPtrAngle angle_; //!< Angular distribution
  UPtrDist energy_; //!< Energy distribution
};


class SourceFile : public SourceDistribution {
public:
  // Constructors
  explicit SourceFile(std::string path);

  // Methods
  Particle::Bank sample(uint64_t* seed) override;

private:
  std::vector<Particle::Bank> sites_; //!< Source sites from a file
};

//==============================================================================
//! Wrapper for custom sources that manages opening/closing shared library
//==============================================================================

class CustomSourceWrapper : public SourceDistribution {
public:
  // Constructors, destructors
  CustomSourceWrapper(std::string path, std::string parameters);
  ~CustomSourceWrapper();

  // Defer implementation to custom source library
  Particle::Bank sample(uint64_t* seed) override
  {
    return custom_source_->sample(seed);
  }

  double strength() const override { return custom_source_->strength(); }
private:
  void* shared_library_; //!< library from dlopen
  std::unique_ptr<SourceDistribution> custom_source_;
};

typedef std::unique_ptr<SourceDistribution> create_custom_source_t(std::string parameters);

//==============================================================================
// Functions
//==============================================================================

//! Initialize source bank from file/distribution
extern "C" void initialize_source();

//! Sample a site from all external source distributions in proportion to their
//! source strength
//! \param[inout] seed Pseudorandom seed pointer
//! \return Sampled source site
Particle::Bank sample_external_source(uint64_t* seed);

void free_memory_source();

} // namespace openmc

#endif // OPENMC_SOURCE_H
