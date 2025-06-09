//! \file source.h
//! \brief External source distributions

#ifndef OPENMC_SOURCE_H
#define OPENMC_SOURCE_H

#include <limits>
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

// Probability distribution for selecting external sources
extern DiscreteIndex external_sources_probability;

} // namespace model

//==============================================================================
//! Abstract source interface
//
//! The Source class provides the interface that must be implemented by derived
//! classes, namely the sample() method that returns a sampled source site. From
//! this base class, source rejection is handled within the
//! sample_with_constraints() method. However, note that some classes directly
//! check for constraints for efficiency reasons (like IndependentSource), in
//! which case the constraints_applied() method indicates that constraints
//! should not be checked a second time from the base class.
//==============================================================================

class Source {
public:
  // Domain types
  enum class DomainType { UNIVERSE, MATERIAL, CELL };
  // Constructors, destructors
  Source() = default;
  explicit Source(pugi::xml_node node);
  virtual ~Source() = default;

  // Methods that can be overridden
  virtual double strength() const { return strength_; }

  //! Sample a source site and apply constraints
  //
  //! \param[inout] seed Pseudorandom seed pointer
  //! \return Sampled site
  SourceSite sample_with_constraints(uint64_t* seed) const;

  //! Sample a source site (without applying constraints)
  //
  //! Sample from the external source distribution
  //! \param[inout] seed Pseudorandom seed pointer
  //! \return Sampled site
  virtual SourceSite sample(uint64_t* seed) const = 0;

  static unique_ptr<Source> create(pugi::xml_node node);

protected:
  // Strategy used for rejecting sites when constraints are applied. KILL means
  // that sites are always accepted but if they don't satisfy constraints, they
  // are given weight 0. RESAMPLE means that a new source site will be sampled
  // until constraints are met.
  enum class RejectionStrategy { KILL, RESAMPLE };

  // Indicates whether derived class already handles constraints
  virtual bool constraints_applied() const { return false; }

  // Methods for constraints
  void read_constraints(pugi::xml_node node);
  bool satisfies_spatial_constraints(Position r) const;
  bool satisfies_energy_constraints(double E) const;
  bool satisfies_time_constraints(double time) const;

  // Data members
  double strength_ {1.0};                  //!< Source strength
  std::unordered_set<int32_t> domain_ids_; //!< Domains to reject from
  DomainType domain_type_;                 //!< Domain type for rejection
  std::pair<double, double> time_bounds_ {-std::numeric_limits<double>::max(),
    std::numeric_limits<double>::max()}; //!< time limits
  std::pair<double, double> energy_bounds_ {
    0, std::numeric_limits<double>::max()}; //!< energy limits
  bool only_fissionable_ {
    false}; //!< Whether site must be in fissionable material
  RejectionStrategy rejection_strategy_ {
    RejectionStrategy::RESAMPLE}; //!< Procedure for rejecting
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

  // Make observing pointers available
  SpatialDistribution* space() const { return space_.get(); }
  UnitSphereDistribution* angle() const { return angle_.get(); }
  Distribution* energy() const { return energy_.get(); }
  Distribution* time() const { return time_.get(); }

  // Make domain type and ids available
  DomainType domain_type() const { return domain_type_; }
  const std::unordered_set<int32_t>& domain_ids() const { return domain_ids_; }

  // Setter for spatial distribution
  void set_space(UPtrSpace space) { space_ = std::move(space); }

protected:
  // Indicates whether derived class already handles constraints
  bool constraints_applied() const override { return true; }

private:
  // Data members
  ParticleType particle_ {ParticleType::neutron}; //!< Type of particle emitted
  UPtrSpace space_;                               //!< Spatial distribution
  UPtrAngle angle_;                               //!< Angular distribution
  UPtrDist energy_;                               //!< Energy distribution
  UPtrDist time_;                                 //!< Time distribution
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
  void load_sites_from_file(
    const std::string& path); //!< Load source sites from file

protected:
  SourceSite sample(uint64_t* seed) const override;

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

  double strength() const override { return compiled_source_->strength(); }

  void setup(const std::string& path, const std::string& parameters);

protected:
  // Defer implementation to custom source library
  SourceSite sample(uint64_t* seed) const override
  {
    return compiled_source_->sample(seed);
  }

private:
  void* shared_library_; //!< library from dlopen
  unique_ptr<Source> compiled_source_;
};

typedef unique_ptr<Source> create_compiled_source_t(std::string parameters);

//==============================================================================
//! Mesh-based source with different distributions for each element
//==============================================================================

// Helper class to sample spatial position on a single mesh element
class MeshElementSpatial : public SpatialDistribution {
public:
  MeshElementSpatial(const Mesh& mesh, int elem_index)
    : mesh_(mesh), elem_index_(elem_index)
  {}

  //! Sample a position from the distribution
  //! \param seed Pseudorandom number seed pointer
  //! \return Sampled position
  Position sample(uint64_t* seed) const override
  {
    return mesh_.sample_element(elem_index_, seed);
  }

private:
  const Mesh& mesh_; //! Reference to mesh
  int elem_index_;   //! Index of mesh element
};

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
  const unique_ptr<IndependentSource>& source(int32_t i) const
  {
    return sources_.size() == 1 ? sources_[0] : sources_[i];
  }

private:
  // Data members
  unique_ptr<MeshSpatial> space_;                 //!< Mesh spatial
  vector<unique_ptr<IndependentSource>> sources_; //!< Source distributions
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
