#ifndef OPENMC_TALLIES_FILTER_H
#define OPENMC_TALLIES_FILTER_H

#include <cstdint>
#include <string>
#include <vector>

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
};

} // namespace openmc

// Without an explicit instantiation of vector<FilterMatch>, the Intel compiler
// will complain about the threadprivate directive on filter_matches. Note that
// this has to happen *outside* of the openmc namespace
template class std::vector<openmc::FilterMatch>;

namespace openmc {

//==============================================================================
//! Modifies tally score events.
//==============================================================================

class Filter
{
public:
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

  virtual void initialize() {}

  int n_bins_;
};

//==============================================================================
// Global variables
//==============================================================================

extern "C" int32_t n_filters;

extern std::vector<FilterMatch> filter_matches;
#pragma omp threadprivate(filter_matches)

extern std::vector<Filter*> tally_filters;

//==============================================================================

extern "C" void free_memory_tally_c();

//==============================================================================

// Filter-related Fortran functions that will be called from C++
extern "C" int verify_filter(int32_t index);
extern "C" Filter* filter_from_f(int32_t index);
extern "C" void filter_update_n_bins(int32_t index);

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_H
