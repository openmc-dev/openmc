#ifndef OPENMC_TALLY_FILTER_H
#define OPENMC_TALLY_FILTER_H

#include <cstdint>
#include <string>
#include <vector>

#include "openmc/hdf5_interface.h"
#include "openmc/particle.h"
#include "openmc/xml_interface.h"


namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

extern "C" int32_t n_filters;

class TallyFilterMatch;
extern std::vector<TallyFilterMatch> filter_matches;
#pragma omp threadprivate(filter_matches)

class TallyFilter;
extern std::vector<TallyFilter*> tally_filters;

//==============================================================================
//! Stores bins and weights for filtered tally events.
//==============================================================================

class TallyFilterMatch
{
public:
  //int i_bin_;
  std::vector<int> bins_;
  std::vector<double> weights_;
  //bool bins_present_;
};

//==============================================================================
//! Modifies tally score events.
//==============================================================================

class TallyFilter
{
public:
  virtual std::string type() const = 0;

  virtual ~TallyFilter() = 0;

  //! Uses an XML input to fill the filter's data fields.
  virtual void from_xml(pugi::xml_node node) = 0;

  //! Matches a tally event to a set of filter bins and weights.
  //!
  //! \param[out] match will contain the matching bins and corresponding
  //!   weights; note that there may be zero matching bins
  virtual void
  get_all_bins(Particle* p, int estimator, TallyFilterMatch& match) const = 0;

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

inline TallyFilter::~TallyFilter() {}

//==============================================================================

extern "C" void free_memory_tally_c();

} // namespace openmc
#endif // OPENMC_TALLY_FILTER_H
