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
  //int i_bin;
  std::vector<int> bins;
  std::vector<double> weights;
  //bool bins_present;
};

//==============================================================================
//! Modifies tally score events.
//==============================================================================

class TallyFilter
{
public:
  virtual std::string type() const = 0;

  virtual ~TallyFilter() = 0;

  virtual void from_xml(pugi::xml_node node) = 0;

  virtual void
  get_all_bins(Particle* p, int estimator, TallyFilterMatch& match) const = 0;

  virtual void
  to_statepoint(hid_t filter_group) const
  {
    write_dataset(filter_group, "type", type());
    write_dataset(filter_group, "n_bins", n_bins_);
  }

  virtual std::string text_label(int bin) const = 0;

  virtual void initialize() {}

  int n_bins_;
};

inline TallyFilter::~TallyFilter() {}

//==============================================================================

extern "C" void free_memory_tally_c();

} // namespace openmc
#endif // OPENMC_TALLY_FILTER_H
