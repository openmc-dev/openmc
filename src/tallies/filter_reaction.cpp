#include "openmc/tallies/filter_reaction.h"

#include <fmt/core.h>

#include "openmc/capi.h"
#include "openmc/endf.h"
#include "openmc/error.h"
#include "openmc/reaction.h"
#include "openmc/xml_interface.h"

namespace openmc {

//==============================================================================
// ReactionFilter implementation
//==============================================================================

void ReactionFilter::from_xml(pugi::xml_node node)
{
  // Read bins as reaction name strings
  auto bins_str = get_node_array<std::string>(node, "bins");

  // Convert reaction names to MT numbers
  vector<int> bins_mt;
  bins_mt.reserve(bins_str.size());
  for (const auto& name : bins_str) {
    bins_mt.push_back(reaction_mt(name));
  }

  this->set_bins(bins_mt);
}

void ReactionFilter::set_bins(span<const int> bins)
{
  // Clear existing bins
  bins_.clear();
  bins_.reserve(bins.size());

  // Copy bins and build lookup map
  for (int64_t i = 0; i < bins.size(); ++i) {
    bins_.push_back(bins[i]);
  }

  n_bins_ = bins_.size();
}

void ReactionFilter::get_all_bins(
  const Particle& p, TallyEstimator estimator, FilterMatch& match) const
{
  // Get the event MT number from the particle
  int event_mt = p.event_mt();

  // Check each bin, considering summation rules
  for (int64_t i = 0; i < bins_.size(); ++i) {
    if (mt_matches(event_mt, bins_[i])) {
      match.bins_.push_back(i);
      match.weights_.push_back(1.0);
    }
  }
}

void ReactionFilter::to_statepoint(hid_t filter_group) const
{
  Filter::to_statepoint(filter_group);

  // Write bins as reaction name strings for human readability
  vector<std::string> names;
  names.reserve(bins_.size());
  for (auto mt : bins_) {
    names.push_back(reaction_name(mt));
  }
  write_dataset(filter_group, "bins", names);
}

std::string ReactionFilter::text_label(int bin) const
{
  return reaction_name(bins_[bin]);
}

} // namespace openmc
