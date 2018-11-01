#ifndef OPENMC_TALLIES_FILTER_DELAYEDGROUP_H
#define OPENMC_TALLIES_FILTER_DELAYEDGROUP_H

#include <vector>

#include "openmc/tallies/filter.h"

namespace openmc {

//==============================================================================
//! Bins outgoing fission neutrons in their delayed groups.
//!
//! The get_all_bins functionality is not actually used.  The bins are manually
//! iterated over in the scoring subroutines.
//==============================================================================

class DelayedGroupFilter : public Filter
{
public:
  ~DelayedGroupFilter() = default;

  std::string type() const override {return "delayedgroup";}

  void from_xml(pugi::xml_node node) override;

  void get_all_bins(const Particle* p, int estimator, FilterMatch& match)
  const override;

  void to_statepoint(hid_t filter_group) const override;

  std::string text_label(int bin) const override;

  std::vector<int> groups_;
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_DELAYEDGROUP_H
