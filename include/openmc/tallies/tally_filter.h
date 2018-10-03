#ifndef OPENMC_TALLY_FILTER_H
#define OPENMC_TALLY_FILTER_H

#include <cstdint>
#include <vector>

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
  int i_bin;
  std::vector<int> bins;
  std::vector<double> weights;
  bool bins_present;
};

//==============================================================================
//! Modifies tally score events.
//==============================================================================

class TallyFilter
{
public:
  virtual void from_xml(pugi::xml_node node) = 0;

  virtual void
  get_all_bins(Particle* p, int estimator, TallyFilterMatch& match) const = 0;

  virtual void initialize() {}
};

//==============================================================================

extern "C" void free_memory_tally_c();

} // namespace openmc
#endif // OPENMC_TALLY_FILTER_H
