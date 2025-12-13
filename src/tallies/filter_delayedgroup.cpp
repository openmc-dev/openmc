#include "openmc/tallies/filter_delayedgroup.h"

#include "openmc/error.h"
#include "openmc/xml_interface.h"

namespace openmc {

void DelayedGroupFilter::from_xml(pugi::xml_node node)
{
  auto groups = get_node_array<int>(node, "bins");
  this->set_groups(groups);
}

void DelayedGroupFilter::set_groups(span<int> groups)
{
  // Clear existing groups
  groups_.clear();
  groups_.reserve(groups.size());

  // Make sure all the group index values are valid.
  // TODO: do these need to be decremented for zero-based indexing?
  for (auto group : groups) {
    if (group < 1) {
      throw std::invalid_argument {"Encountered delayedgroup bin with index " +
                                   std::to_string(group) +
                                   " which is less than 1"};
    } else if (group > MAX_DELAYED_GROUPS) {
      throw std::invalid_argument {
        "Encountered delayedgroup bin with index " + std::to_string(group) +
        " which is greater than MAX_DELAYED_GROUPS (" +
        std::to_string(MAX_DELAYED_GROUPS) + ")"};
    }
    groups_.push_back(group);
  }

  n_bins_ = groups_.size();
}

void DelayedGroupFilter::get_all_bins(
  const Particle& p, TallyEstimator estimator, FilterMatch& match) const
{
  // Note that the bin is set to zero here, but bins outside zero are
  // tallied to regardless. This is because that logic has to be handled
  // in the scoring code instead where looping over the delayed
  // group takes place (tally_scoring.cpp).
  match.bins_.push_back(0);
  match.weights_.push_back(1.0);
}

void DelayedGroupFilter::to_statepoint(hid_t filter_group) const
{
  Filter::to_statepoint(filter_group);
  write_dataset(filter_group, "bins", groups_);
}

std::string DelayedGroupFilter::text_label(int bin) const
{
  return "Delayed Group " + std::to_string(groups_[bin]);
}

} // namespace openmc
