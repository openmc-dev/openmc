#ifndef OPENMC_TALLIES_FILTER_DELAYEDGROUP_H
#define OPENMC_TALLIES_FILTER_DELAYEDGROUP_H

#include <vector>

#include <gsl/gsl>

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
  //----------------------------------------------------------------------------
  // Constructors, destructors

  ~DelayedGroupFilter() = default;

  //----------------------------------------------------------------------------
  // Methods

  std::string type() const override {return "delayedgroup";}

  void from_xml(pugi::xml_node node) override;

  void get_all_bins(const Particle* p, TallyEstimator estimator, FilterMatch& match)
  const override;

  void to_statepoint(hid_t filter_group) const override;

  std::string text_label(int bin) const override;

  //----------------------------------------------------------------------------
  // Accessors

  const std::vector<int>& groups() const { return groups_; }

  void set_groups(gsl::span<int> groups);

private:
  //----------------------------------------------------------------------------
  // Data members

  std::vector<int> groups_;
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_DELAYEDGROUP_H
