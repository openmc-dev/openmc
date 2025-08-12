#ifndef OPENMC_TALLIES_TALLY_H
#define OPENMC_TALLIES_TALLY_H

#include "openmc/constants.h"
#include "openmc/memory.h" // for unique_ptr
#include "openmc/span.h"
#include "openmc/tallies/filter.h"
#include "openmc/tallies/trigger.h"
#include "openmc/vector.h"

#include "pugixml.hpp"
#include "xtensor/xfixed.hpp"
#include "xtensor/xtensor.hpp"

#include <string>
#include <unordered_map>

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

  int id() const { return id_; }

  void set_active(bool active) { active_ = active; }

  void set_multiply_density(bool value) { multiply_density_ = value; }

  void set_writable(bool writable) { writable_ = writable; }

  void set_scores(pugi::xml_node node);

  void set_scores(const vector<std::string>& scores);

  std::vector<std::string> scores() const;

  int32_t n_scores() const { return scores_.size(); }

  void set_nuclides(pugi::xml_node node);

  void set_nuclides(const vector<std::string>& nuclides);

  const xt::xtensor<double, 3>& results() const { return results_; }

  //! returns vector of indices corresponding to the tally this is called on
  const vector<int32_t>& filters() const { return filters_; }

  //! returns a vector of filter types for the tally
  std::vector<FilterType> filter_types() const;

  //! returns a mapping of filter types to index into the tally's filters
  std::unordered_map<FilterType, int32_t> filter_indices() const;

  //! \brief Returns the tally filter at index i
  int32_t filters(int i) const { return filters_[i]; }

  //! \brief Return a const pointer to a filter instance based on type. Always
  //! returns the first matching filter type
  template<class T>
  const T* get_filter() const
  {
    const T* out;
    for (auto filter_idx : filters_) {
      if ((out = dynamic_cast<T*>(model::tally_filters[filter_idx].get())))
        return out;
    }
    return nullptr;
  }

  template<class T>
  const T* get_filter(int idx) const
  {
    if (const T* out = dynamic_cast<T*>(model::tally_filters[filters_.at(idx)]))
      return out;
    return nullptr;
  }

  //! \brief Check if this tally has a specified type of filter
  bool has_filter(FilterType filter_type) const;

  void set_filters(span<Filter*> filters);

  //! Given already-set filters, set the stride lengths
  void set_strides();

  int32_t strides(int i) const { return strides_[i]; }

  int32_t n_filter_bins() const { return n_filter_bins_; }

  bool multiply_density() const { return multiply_density_; }

  bool writable() const { return writable_; }

  bool vov_results() const { return vov_; }

  //----------------------------------------------------------------------------
  // Other methods.

  void add_filter(Filter* filter);

  void init_triggers(pugi::xml_node node);

  void init_results();

  void reset();

  void accumulate();

  void accumulate_vov();

  //! return the index of a score specified by name
  int score_index(const std::string& score) const;

  //! Tally results reshaped according to filter sizes
  xt::xarray<double> get_reshaped_data() const;

  //! A string representing the i-th score on this tally
  std::string score_name(int score_idx) const;

  //! A string representing the i-th nuclide on this tally
  std::string nuclide_name(int nuclide_idx) const;

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
  xt::xtensor<double, 3> results_;

  //! True if this tally should be written to statepoint files
  bool writable_ {true};

  //----------------------------------------------------------------------------
  // Miscellaneous public members.

  // We need to have quick access to some filters.  The following gives indices
  // for various filters that could be in the tally or C_NONE if they are not
  // present.
  int energyout_filter_ {C_NONE};
  int delayedgroup_filter_ {C_NONE};

  vector<Trigger> triggers_;

  int deriv_ {C_NONE}; //!< Index of a TallyDerivative object for diff tallies.

private:
  //----------------------------------------------------------------------------
  // Private data.

  vector<int32_t> filters_; //!< Filter indices in global filters array

  //! Index strides assigned to each filter to support 1D indexing.
  vector<int32_t> strides_;

  int32_t n_filter_bins_ {0};

  //! Whether to multiply by atom density for reaction rates
  bool multiply_density_ {true};

  bool vov_ = {false};

  int64_t index_;
};

//==============================================================================
// Global variable declarations
//==============================================================================

namespace model {
extern std::unordered_map<int, int> tally_map;
extern vector<unique_ptr<Tally>> tallies;
extern vector<int> active_tallies;
extern vector<int> active_analog_tallies;
extern vector<int> active_tracklength_tallies;
extern vector<int> active_collision_tallies;
extern vector<int> active_meshsurf_tallies;
extern vector<int> active_surface_tallies;
extern vector<int> active_pulse_height_tallies;
extern vector<int> pulse_height_cells;
} // namespace model

namespace simulation {
//! Global tallies (such as k-effective estimators)
extern xt::xtensor_fixed<double, xt::xshape<N_GLOBAL_TALLIES, 3>>
  global_tallies;

//! Number of realizations for global tallies
extern "C" int32_t n_realizations;
} // namespace simulation

extern double global_tally_absorption;
extern double global_tally_collision;
extern double global_tally_tracklength;
extern double global_tally_leakage;

//==============================================================================
// Non-member functions
//==============================================================================

//! Read tally specification from tallies.xml
void read_tallies_xml();

//! Read tally specification from an XML node
//! \param[in] root node of tallies XML element
void read_tallies_xml(pugi::xml_node root);

//! \brief Accumulate the sum of the contributions from each history within the
//! batch to a new random variable
void accumulate_tallies();

//! Determine which tallies should be active
void setup_active_tallies();

// Alias for the type returned by xt::adapt(...). N is the dimension of the
// multidimensional array
template<std::size_t N>
using adaptor_type =
  xt::xtensor_adaptor<xt::xbuffer_adaptor<double*&, xt::no_ownership>, N>;

#ifdef OPENMC_MPI
//! Collect all tally results onto master process
void reduce_tally_results();
#endif

void free_memory_tally();

} // namespace openmc

#endif // OPENMC_TALLIES_TALLY_H
