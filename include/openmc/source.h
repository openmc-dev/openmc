//! \file source.h
//! \brief External source distributions

#ifndef OPENMC_SOURCE_H
#define OPENMC_SOURCE_H

#include <string>

#include "pugixml.hpp"

#include "openmc/distribution_multi.h"
#include "openmc/distribution_spatial.h"
#include "openmc/memory.h"
#include "openmc/particle.h"
#include "openmc/vector.h"

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

class Source;

namespace model {

extern vector<unique_ptr<Source>> external_sources;

} // namespace model

#ifdef __CUDACC__
namespace gpu {
extern __constant__ unique_ptr<Source>* external_sources;
extern __constant__ int n_external_sources;
} // namespace gpu
#endif

//==============================================================================
//! Abstract source interface
//==============================================================================

class Source {
public:
  virtual ~Source() = default;
  Source() = default;
  Source(Source&&) = default;

  // Methods that must be implemented. The particle pointer, if not null,
  // can be used to not need to allocate a particle object as is required
  // by openmc_find_cell, which can be time-consuming in the rejection loop.
  HD virtual SourceSite sample(uint64_t* seed, Particle* p = nullptr) const = 0;

  // Methods that can be overridden
  HD virtual double strength() const { return 1.0; }
};

//==============================================================================
//! Source composed of independent spatial, angle, and energy distributions
//==============================================================================

class IndependentSource : public Source {
public:
  // Constructors
  IndependentSource(
    UPtrSpace space, UPtrAngle angle, unique_ptr<Distribution> energy);
  explicit IndependentSource(pugi::xml_node node);
  IndependentSource(IndependentSource&&) = default;

  //! Sample from the external source distribution
  //! \param[inout] seed Pseudorandom seed pointer
  //! \return Sampled site
  HD SourceSite sample(uint64_t* seed, Particle* p) const override;

  // Properties
  HD ParticleType particle_type() const { return particle_; }
  HD double strength() const override { return strength_; }

  // Make observing pointers available
  HD SpatialDistribution* space() const { return space_.get(); }
  HD UnitSphereDistribution* angle() const { return angle_.get(); }
  HD Distribution* energy() const { return energy_.get(); }

private:
  ParticleType particle_ {ParticleType::neutron}; //!< Type of particle emitted
  double strength_ {1.0}; //!< Source strength
  UPtrSpace space_; //!< Spatial distribution
  UPtrAngle angle_; //!< Angular distribution
  unique_ptr<Distribution> energy_; //!< Energy distribution
};

//==============================================================================
//! Source composed of particles read from a file
//==============================================================================

class FileSource : public Source {
public:
  // Constructors
  explicit FileSource(std::string const& path);
  FileSource(FileSource&&) = default;

  // Methods
  HD SourceSite sample(uint64_t* seed, Particle* p) const override;

private:
  vector<SourceSite> sites_; //!< Source sites from a file
};

//==============================================================================
//! Wrapper for custom sources that manages opening/closing shared library
//!
//! This does not work in GPU mode, since the CUDA compiler fails to figure
//! out the requisite stack size for this. Moreover, dynamically  linked
//! device code like this is probably impossible.
//==============================================================================

#ifndef __CUDACC__
class CustomSourceWrapper : public Source {
public:
  // Constructors, destructors
  CustomSourceWrapper(std::string const& path, std::string const& parameters);
  CustomSourceWrapper(CustomSourceWrapper&&) = default;
  virtual ~CustomSourceWrapper();

  // Defer implementation to custom source library
  HD SourceSite sample(uint64_t* seed, Particle* p = nullptr) const override
  {
    return custom_source_->sample(seed);
  }

  HD double strength() const override { return custom_source_->strength(); }

private:
  void* shared_library_; //!< library from dlopen
  unique_ptr<Source> custom_source_;
};

typedef unique_ptr<Source> create_custom_source_t(std::string parameters);
#endif

//==============================================================================
// Functions
//==============================================================================

//! Initialize source bank from file/distribution
extern "C" void initialize_source();

//! Sample a site from all external source distributions in proportion to their
//! source strength
//! \param[inout] seed Pseudorandom seed pointer
//! \return Sampled source site
HD SourceSite sample_external_source(uint64_t* seed, Particle* p = nullptr);

void free_memory_source();

} // namespace openmc

#endif // OPENMC_SOURCE_H
