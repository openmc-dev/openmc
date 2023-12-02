//! \file source.h
//! \brief External source distributions

#ifndef OPENMC_SOURCE_H
#define OPENMC_SOURCE_H

#include <unordered_set>

#include "pugixml.hpp"

#include "openmc/distribution_multi.h"
#include "openmc/distribution_spatial.h"
#include "openmc/memory.h"
#include "openmc/particle.h"
#include "openmc/vector.h"

namespace openmc {

//==============================================================================
// Constants
//==============================================================================

// Maximum number of external source spatial resamples to encounter before an
// error is thrown.
constexpr int EXTSRC_REJECT_THRESHOLD {10000};
constexpr double EXTSRC_REJECT_FRACTION {0.05};

//==============================================================================
// Global variables
//==============================================================================

class Source;

namespace model {

extern vector<unique_ptr<Source>> external_sources;

} // namespace model

//==============================================================================
//! Abstract source interface
//==============================================================================

class Source {
public:
  virtual ~Source() = default;

  // Methods that must be implemented
  virtual SourceSite sample(uint64_t* seed) const = 0;

  // Methods that can be overridden
  virtual double strength() const { return 1.0; }

  static unique_ptr<Source> create(pugi::xml_node node);
};

//==============================================================================
//! Source composed of independent spatial, angle, energy, and time
//! distributions
//==============================================================================

class IndependentSource : public Source {
public:
  // Constructors
  IndependentSource(
    UPtrSpace space, UPtrAngle angle, UPtrDist energy, UPtrDist time);
  explicit IndependentSource(pugi::xml_node node);

  //! Sample from the external source distribution
  //! \param[inout] seed Pseudorandom seed pointer
  //! \return Sampled site
  SourceSite sample(uint64_t* seed) const override;

  // Properties
  ParticleType particle_type() const { return particle_; }
  double strength() const override { return strength_; }

  // Make observing pointers available
  SpatialDistribution* space() const { return space_.get(); }
  UnitSphereDistribution* angle() const { return angle_.get(); }
  Distribution* energy() const { return energy_.get(); }
  Distribution* time() const { return time_.get(); }

private:
  // Domain types
  enum class DomainType { UNIVERSE, MATERIAL, CELL };

  // Data members
  ParticleType particle_ {ParticleType::neutron}; //!< Type of particle emitted
  double strength_ {1.0};                         //!< Source strength
  UPtrSpace space_;                               //!< Spatial distribution
  UPtrAngle angle_;                               //!< Angular distribution
  UPtrDist energy_;                               //!< Energy distribution
  UPtrDist time_;                                 //!< Time distribution
  DomainType domain_type_;                        //!< Domain type for rejection
  std::unordered_set<int32_t> domain_ids_;        //!< Domains to reject from
};

//==============================================================================
//! Source composed of particles read from a file
//==============================================================================

class FileSource : public Source {
public:
  // Constructors
  explicit FileSource(pugi::xml_node node);
  explicit FileSource(const std::string& path);

  // Methods
  SourceSite sample(uint64_t* seed) const override;
  void load_sites_from_file(
    const std::string& path); //!< Load source sites from file
private:
  vector<SourceSite> sites_; //!< Source sites from a file
};

//==============================================================================
//! Wrapper for custom sources that manages opening/closing shared library
//==============================================================================

class CompiledSourceWrapper : public Source {
public:
  // Constructors, destructors
  CompiledSourceWrapper(pugi::xml_node node);
  ~CompiledSourceWrapper();

  // Defer implementation to custom source library
  SourceSite sample(uint64_t* seed) const override
  {
    return compiled_source_->sample(seed);
  }

  double strength() const override { return compiled_source_->strength(); }

  void setup(const std::string& path, const std::string& parameters);

private:
  void* shared_library_; //!< library from dlopen
  unique_ptr<Source> compiled_source_;
};

typedef unique_ptr<Source> create_compiled_source_t(std::string parameters);

//==============================================================================
//! Mesh-based source with different distributions for each element
//==============================================================================

class MeshSource : public Source {
public:
  // Constructors
  explicit MeshSource(pugi::xml_node node);

  //! Sample from the external source distribution
  //! \param[inout] seed Pseudorandom seed pointer
  //! \return Sampled site
  SourceSite sample(uint64_t* seed) const override;

  // Properties
  double strength() const override { return space_->total_strength(); }

  // Accessors
  const std::unique_ptr<Source>& source(int32_t i) const
  {
    return sources_.size() == 1 ? sources_[0] : sources_[i];
  }

private:
  // Data members
  unique_ptr<MeshSpatial> space_;           //!< Mesh spatial
  vector<std::unique_ptr<Source>> sources_; //!< Source distributions
};

//==============================================================================
// Functions
//==============================================================================

//! Initialize source bank from file/distribution
extern "C" void initialize_source();

//! Sample a site from all external source distributions in proportion to their
//! source strength
//! \param[inout] seed Pseudorandom seed pointer
//! \return Sampled source site
SourceSite sample_external_source(uint64_t* seed);

void free_memory_source();

} // namespace openmc

#endif // OPENMC_SOURCE_H
