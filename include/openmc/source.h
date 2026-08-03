//! \file source.h
//! \brief External source distributions

#ifndef OPENMC_SOURCE_H
#define OPENMC_SOURCE_H

#include <atomic>
#include <limits>
#include <unordered_set>
#include <utility> // for pair

#include "pugixml.hpp"

#include "openmc/array.h"
#include "openmc/distribution.h"
#include "openmc/distribution_multi.h"
#include "openmc/distribution_spatial.h"
#include "openmc/memory.h"
#include "openmc/particle_type.h"
#include "openmc/vector.h"

namespace openmc {

//==============================================================================
// Constants
//==============================================================================

// Minimum number of external source sites rejected before checking againts the
// source_rejection_fraction
constexpr int EXTSRC_REJECT_THRESHOLD {10000};

// Maximum number of source rejections allowed while sampling a single site
constexpr int64_t MAX_SOURCE_REJECTIONS_PER_SAMPLE {1'000'000};

//==============================================================================
// Global variables
//==============================================================================

// Cumulative counters for source rejection diagnostics. These are atomic to
// allow thread-safe concurrent sampling of external sources.
extern std::atomic<int64_t> source_n_accept;
extern std::atomic<int64_t> source_n_reject;

class Source;

namespace model {

extern vector<unique_ptr<Source>> external_sources;
extern vector<unique_ptr<Source>> adjoint_sources;

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
  ParticleType particle_; //!< Type of particle emitted
  UPtrSpace space_;       //!< Spatial distribution
  UPtrAngle angle_;       //!< Angular distribution
  UPtrDist energy_;       //!< Energy distribution
  UPtrDist time_;         //!< Time distribution
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
  vector<SourceSite> sites_; //!< Source sites
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
  MeshElementSpatial(int32_t mesh_index, int elem_index)
    : mesh_index_(mesh_index), elem_index_(elem_index)
  {}

  //! Sample a position from the distribution
  //! \param seed Pseudorandom number seed pointer
  //! \return (sampled position, importance weight)
  std::pair<Position, double> sample(uint64_t* seed) const override;

private:
  int32_t mesh_index_ {C_NONE}; //!< Index in global meshes array
  int elem_index_;              //! Index of mesh element
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
//! Parametric tokamak plasma neutron source
//!
//! This source samples neutron positions from a tokamak plasma geometry using
//! Miller-style flux surface parameterization with user-specified emission
//! profiles and energy distributions.
//!
//! Flux surface parameterization:
//!   R = R0 + r*cos(alpha + delta*sin(alpha)) + Delta*(1 - (r/a)^2)
//!   Z = kappa * r * sin(alpha)
//!
//! The sampling algorithm:
//! 1. Sample minor radius r from precomputed CDF of S(r) * Jacobian
//! 2. Sample poloidal angle alpha from conditional P(alpha|r) using mixture
//!    of precomputed CDFs weighted by functions of r
//! 3. Sample energy and time from user-provided distribution(s)
//! 4. Sample isotropic direction
//! 5. Sample toroidal angle phi uniformly in [phi_start, phi_start +
//! phi_extent]
//! 6. Transform (r, alpha, phi) to Cartesian (x, y, z), applying the optional
//!    vertical shift of the plasma center
//!
//! The user provides the emission density S(r) directly (e.g., from transport
//! codes like TRANSP, ASTRA, etc.), allowing full flexibility in reaction
//! physics calculations. S(r) is a profile in arbitrary units sampled on the
//! r_over_a grid; only its shape matters, since it is normalized internally.
//! Energy distributions can be specified as either a single distribution for
//! all r, or one distribution per radial point.
//==============================================================================

class TokamakSource : public Source {
public:
  // Constructors
  explicit TokamakSource(pugi::xml_node node);

  //! Sample from the tokamak source distribution
  //! \param[inout] seed Pseudorandom seed pointer
  //! \return Sampled site
  SourceSite sample(uint64_t* seed) const override;

private:
  //==========================================================================
  // Private methods

  //! Precompute data structures for efficient sampling
  void precompute_sampling_distributions();

  //! Sample minor radius from marginal CDF
  //! \param seed Pseudorandom seed pointer
  //! \return Sampled r/a value
  double sample_r_over_a(uint64_t* seed) const;

  //! Sample poloidal angle given r using mixture of precomputed CDFs
  //! \param r_norm Normalized minor radius r/a
  //! \param seed Pseudorandom seed pointer
  //! \return Sampled poloidal angle alpha [rad]
  double sample_poloidal_angle(double r_norm, uint64_t* seed) const;

  //! Compute the k-th mixture weight w_k(r) * I_hat_k for poloidal sampling
  //! \param k Basis function index (0-5)
  //! \param r Normalized minor radius r/a
  //! \return Mixture weight for component k
  double mixture_weight(int k, double r) const;

  //! Sample energy from the distribution(s)
  //! \param r_norm Normalized minor radius r/a (for distribution selection)
  //! \param seed Pseudorandom seed pointer
  //! \return (Sampled energy [eV], importance weight)
  std::pair<double, double> sample_energy(double r_norm, uint64_t* seed) const;

  //! Transform from flux coordinates (r, alpha, phi) to Cartesian (x, y, z)
  //! \param r Minor radius [cm]
  //! \param alpha Poloidal angle [rad]
  //! \param phi Toroidal angle [rad]
  //! \return Position in Cartesian coordinates [cm]
  Position flux_to_cartesian(double r, double alpha, double phi) const;

  //==========================================================================
  // Data members

  // Emission profile (input)
  vector<double> r_over_a_;         //!< Normalized minor radius grid points
  vector<double> emission_density_; //!< Emission density S(r) at grid points

  // Energy distribution(s): either 1 for all r, or one per r point
  vector<unique_ptr<Distribution>> energy_dists_;

  // Time distribution (defaults to a delta distribution at t=0)
  UPtrDist time_;

  // Angular distribution (isotropic)
  UPtrAngle angle_;

  // Tokamak geometry parameters
  double major_radius_;    //!< Major radius R0 [cm]
  double minor_radius_;    //!< Minor radius a [cm]
  double elongation_;      //!< Elongation kappa
  double triangularity_;   //!< Triangularity delta
  double shafranov_shift_; //!< Shafranov shift Delta [cm]
  double vertical_shift_;  //!< Vertical shift of plasma center [cm]

  // Normalized geometry parameters (precomputed for efficiency)
  double epsilon_;     //!< Inverse aspect ratio a/R0
  double delta_tilde_; //!< Normalized Shafranov shift Delta/a

  // Toroidal angle bounds
  double phi_start_;  //!< Starting toroidal angle [rad]
  double phi_extent_; //!< Toroidal angle extent [rad]

  // Precomputed distribution for radial sampling
  unique_ptr<Tabular> radial_dist_;

  // Coefficients of the radial geometric polynomial: A*r - B*r^2 - C*r^3
  // Also used as the analytical normalization for poloidal mixture weights
  double radial_poly_a_; //!< 1 + ε·Δ̃
  double radial_poly_b_; //!< (3/8)·c₁·ε
  double radial_poly_c_; //!< 2·ε·Δ̃

  // Precomputed Bernstein basis functions for poloidal angle sampling.
  // Using the factorization f(r_tilde, alpha) = R_tilde x J_tilde where:
  //   R_tilde = b0*(1-r)^2 + 2*b1*r*(1-r) + b2*r^2  (Bernstein quadratic)
  //   J_tilde = b3*(1-r) + b4*r                      (Bernstein linear)
  // The product gives 6 non-negative basis functions g_k(alpha) with
  // weights w_k(r_tilde) that are products of Bernstein polynomials.
  // Distributions are tabulated on [0, pi] exploiting up-down symmetry.
  static constexpr int N_POLOIDAL_BASIS = 6; //!< Number of basis functions
  int n_alpha_; //!< Number of poloidal angle grid points
  array<unique_ptr<Tabular>, N_POLOIDAL_BASIS>
    poloidal_dists_; //!< Distributions for each basis function g_k(alpha)
  array<double, N_POLOIDAL_BASIS>
    poloidal_integrals_; //!< Integrals of g_k(alpha) over [0, pi]
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

//! Reset cumulative source rejection counters
void reset_source_rejection_counters();

} // namespace openmc

#endif // OPENMC_SOURCE_H
