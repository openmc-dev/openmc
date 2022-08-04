#ifndef OPENMC_TALLIES_TALLY_H
#define OPENMC_TALLIES_TALLY_H

#include "openmc/constants.h"
#include "openmc/tallies/filter.h"
#include "openmc/tallies/trigger.h"
#include "openmc/vector.h"

#include <gsl/gsl>
#include "pugixml.hpp"
#include "xtensor/xfixed.hpp"
#include "xtensor/xtensor.hpp"

#include <array>
#include <memory> // for unique_ptr
#include <unordered_map>
#include <string>
#include <vector>

namespace openmc {

//==============================================================================
//! A user-specified flux-weighted (or current) measurement.
//==============================================================================

class Tally {
public:
  //----------------------------------------------------------------------------
  // Constructors, destructors, factory functions
  explicit Tally(int32_t id);
  explicit Tally(pugi::xml_node node);
  ~Tally();
  static Tally* create(int32_t id = -1);

  //----------------------------------------------------------------------------
  // Accessors

  void set_id(int32_t id);

  void set_active(bool active) { active_ = active; }

  void set_writable(bool writable) { writable_ = writable; }

  void set_scores(pugi::xml_node node);

  void set_scores(const std::vector<std::string>& scores);

  void set_nuclides(pugi::xml_node node);

  void set_nuclides(const std::vector<std::string>& nuclides);

  const std::vector<int32_t>& filters() const {return filters_;}

  int32_t filters(int i) const {return filters_[i];}

  void set_filters(gsl::span<Filter*> filters);

  int32_t strides(int i) const {return strides_[i];}

  int32_t n_filter_bins() const {return n_filter_bins_;}

  bool writable() const { return writable_;}

  //----------------------------------------------------------------------------
  // Other methods.

  void add_filter(Filter* filter) { set_filters({&filter, 1}); }

  void init_triggers(pugi::xml_node node);

  void init_results();

  void reset();

  void accumulate();

  //! A string representing the i-th score on this tally
  std::string score_name(int score_idx) const;

  void copy_to_device();
  void release_from_device();

  //----------------------------------------------------------------------------
  // Major public data members.

  int id_ {C_NONE}; //!< User-defined identifier

  std::string name_; //!< User-defined name

  TallyType type_ {TallyType::VOLUME}; //!< e.g. volume, surface current

  //! Event type that contributes to this tally
  TallyEstimator estimator_ {TallyEstimator::TRACKLENGTH};

  //! Whether this tally is currently being updated
  bool active_ {false};

  //! Number of realizations
  int n_realizations_ {0};

  vector<int> scores_; //!< Filter integrands (e.g. flux, fission)

  //! Index of each nuclide to be tallied.  -1 indicates total material.
  vector<int> nuclides_ {-1};

  //! Results for each bin -- the first dimension of the array is for the
  //! combination of filters (e.g. specific cell, specific energy group, etc.)
  //! and the second dimension of the array is for scores (e.g. flux, total
  //! reaction rate, fission reaction rate, etc.)
  double* results_;
  size_t results_size_ {0};
  size_t n_scores_;

  #pragma omp declare target
  const double* results(gsl::index i, gsl::index j, TallyResult k) const;
  double* results(gsl::index i, gsl::index j, TallyResult k);
  std::array<size_t, 3> results_shape() const;
  #pragma omp end declare target

  //! True if this tally should be written to statepoint files
  bool writable_ {true};

  //----------------------------------------------------------------------------
  // Miscellaneous public members.

  // We need to have quick access to some filters.  The following gives indices
  // for various filters that could be in the tally or C_NONE if they are not
  // present.
  int energyout_filter_ {C_NONE};
  int delayedgroup_filter_ {C_NONE};

  std::vector<Trigger> triggers_;

  int deriv_ {C_NONE}; //!< Index of a TallyDerivative object for diff tallies.

private:
  //----------------------------------------------------------------------------
  // Private data.

  std::vector<int32_t> filters_; //!< Filter indices in global filters array

  //! Index strides assigned to each filter to support 1D indexing.
  std::vector<int32_t> strides_;

  int32_t n_filter_bins_ {0};

  gsl::index index_;
};

//==============================================================================
// Global variable declarations
//==============================================================================

namespace model {
  extern std::unordered_map<int, int> tally_map;
  extern std::vector<int> active_tallies;
  extern std::vector<int> active_analog_tallies;
  extern std::vector<int> active_tracklength_tallies;
  extern std::vector<int> active_collision_tallies;
  extern std::vector<int> active_meshsurf_tallies;
  extern std::vector<int> active_surface_tallies;

  #pragma omp declare target
  extern Tally* tallies;
  extern size_t tallies_size;
  extern int* device_active_tallies;
  extern size_t active_tallies_size;
  extern int* device_active_collision_tallies;
  extern size_t active_collision_tallies_size;
  extern int* device_active_tracklength_tallies;
  extern size_t active_tracklength_tallies_size;
  #pragma omp end declare target
}

namespace simulation {
  //! Global tallies (such as k-effective estimators)
  extern xt::xtensor_fixed<double, xt::xshape<N_GLOBAL_TALLIES, 3>> global_tallies;

  //! Number of realizations for global tallies
  extern "C" int32_t n_realizations;
}

#pragma omp declare target
extern double global_tally_absorption;
extern double global_tally_collision;
extern double global_tally_tracklength;
extern double global_tally_leakage;
#pragma omp end declare target

//==============================================================================
// Non-member functions
//==============================================================================

//! Read tally specification from tallies.xml
void read_tallies_xml();

//! \brief Accumulate the sum of the contributions from each history within the
//! batch to a new random variable
void accumulate_tallies();

//! Determine which tallies should be active
void setup_active_tallies();

// Alias for the type returned by xt::adapt(...). N is the dimension of the
// multidimensional array
template <std::size_t N>
using adaptor_type = xt::xtensor_adaptor<xt::xbuffer_adaptor<double*&, xt::no_ownership>, N>;

#ifdef OPENMC_MPI
//! Collect all tally results onto master process
void reduce_tally_results();
#endif

void free_memory_tally();

} // namespace openmc

#endif // OPENMC_TALLIES_TALLY_H
