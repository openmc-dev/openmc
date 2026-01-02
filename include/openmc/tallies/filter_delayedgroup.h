#ifndef OPENMC_TALLIES_FILTER_DELAYEDGROUP_H
#define OPENMC_TALLIES_FILTER_DELAYEDGROUP_H

#include "openmc/span.h"
#include "openmc/tallies/filter.h"
#include "openmc/vector.h"

namespace openmc {

//==============================================================================
//! Bins outgoing fission neutrons in their delayed groups.
//!
//! The get_all_bins functionality is not actually used.  The bins are manually
//! iterated over in the scoring subroutines.
//==============================================================================

class DelayedGroupFilter : public Filter {
public:
  //----------------------------------------------------------------------------
  // Constructors, destructors

  ~DelayedGroupFilter() = default;

  //----------------------------------------------------------------------------
  // Methods

  std::string type_str() const override { return "delayedgroup"; }
  FilterType type() const override { return FilterType::DELAYED_GROUP; }

  void from_xml(pugi::xml_node node) override;

  void get_all_bins(const Particle& p, TallyEstimator estimator,
    FilterMatch& match) const override;

  void to_statepoint(hid_t filter_group) const override;

  std::string text_label(int bin) const override;

  //----------------------------------------------------------------------------
  // Accessors

  const vector<int>& groups() const { return groups_; }

  void set_groups(span<int> groups);

private:
  //----------------------------------------------------------------------------
  // Data members

  vector<int> groups_;
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_DELAYEDGROUP_H
