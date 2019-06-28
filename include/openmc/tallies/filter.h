#ifndef OPENMC_TALLIES_FILTER_H
#define OPENMC_TALLIES_FILTER_H

#include <cstdint>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include <gsl/gsl>

#include "openmc/hdf5_interface.h"
#include "openmc/particle.h"
#include "pugixml.hpp"


namespace openmc {

//==============================================================================
//! Stores bins and weights for filtered tally events.
//==============================================================================

class FilterMatch
{
public:
  std::vector<int> bins_;
  std::vector<double> weights_;
  int i_bin_;
  bool bins_present_ {false};
};

} // namespace openmc

// Without an explicit instantiation of vector<FilterMatch>, the Intel compiler
// will complain about the threadprivate directive on filter_matches. Note that
// this has to happen *outside* of the openmc namespace
extern template class std::vector<openmc::FilterMatch>;

namespace openmc {

//==============================================================================
//! Modifies tally score events.
//==============================================================================

class Filter
{
public:
  // Default constructor
  Filter();

  //! Create a new tally filter
  //
  //! \param[in] type  Type of the filter
  //! \param[in] id  Unique ID for the filter. If none is passed, an ID is
  //!    automatically assigned
  //! \return Pointer to the new filter object
  static Filter* create(const std::string& type, int32_t id = -1);

  //! Create a new tally filter from an XML node
  //
  //! \param[in] node XML node
  //! \return Pointer to the new filter object
  static Filter* create(pugi::xml_node node);

  virtual ~Filter() = default;

  virtual std::string type() const = 0;

  //! Uses an XML input to fill the filter's data fields.
  virtual void from_xml(pugi::xml_node node) = 0;

  //! Matches a tally event to a set of filter bins and weights.
  //!
  //! \param[out] match will contain the matching bins and corresponding
  //!   weights; note that there may be zero matching bins
  virtual void
  get_all_bins(const Particle* p, int estimator, FilterMatch& match) const = 0;

  //! Assign a unique ID to the filter
  //
  //! \param[in]  Unique ID to assign
  void set_id(int32_t id);

  gsl::index index() const { return index_; }

  //! Writes data describing this filter to an HDF5 statepoint group.
  virtual void
  to_statepoint(hid_t filter_group) const
  {
    write_dataset(filter_group, "type", type());
    write_dataset(filter_group, "n_bins", n_bins_);
  }

  //! Return a string describing a filter bin for the tallies.out file.
  //
  //! For example, an `EnergyFilter` might return the string
  //! "Incoming Energy [0.625E-6, 20.0)".
  virtual std::string text_label(int bin) const = 0;

  int32_t id_ {-1};

  int n_bins_;
private:
  gsl::index index_;
};

//==============================================================================
// Global variables
//==============================================================================

namespace simulation {

extern std::vector<FilterMatch> filter_matches;
#pragma omp threadprivate(filter_matches)

} // namespace simulation

namespace model {
  extern "C" int32_t n_filters;
  extern std::vector<std::unique_ptr<Filter>> tally_filters;
  extern std::unordered_map<int, int> filter_map;
}

//==============================================================================
// Non-member functions
//==============================================================================

//! Make sure index corresponds to a valid filter
int verify_filter(int32_t index);

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_H
