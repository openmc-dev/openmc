#include "openmc/tallies/filter_delayedgroup.h"

#include "openmc/error.h"
#include "openmc/xml_interface.h"

namespace openmc {

void
DelayedGroupFilter::from_xml(pugi::xml_node node)
{
  groups_ = get_node_array<int>(node, "bins");
  n_bins_ = groups_.size();

  // Make sure all the group index values are valid.
  // TODO: do these need to be decremented for zero-based indexing?
  for (auto group : groups_) {
    if (group < 1) {
      fatal_error("Encountered delayedgroup bin with index "
        + std::to_string(group) + " which is less than 1");
    } else if (group > MAX_DELAYED_GROUPS) {
      fatal_error("Encountered delayedgroup bin with index "
        + std::to_string(group) + " which is greater than MAX_DELATED_GROUPS ("
        + std::to_string(MAX_DELAYED_GROUPS) + ")");
    }
  }
}

void
DelayedGroupFilter::get_all_bins(const Particle* p, int estimator,
                                 FilterMatch& match) const
{
  //TODO: off-by-one
  match.bins_.push_back(1);
  match.weights_.push_back(1.0);
}

void
DelayedGroupFilter::to_statepoint(hid_t filter_group) const
{
  Filter::to_statepoint(filter_group);
  write_dataset(filter_group, "bins", groups_);
}

std::string
DelayedGroupFilter::text_label(int bin) const
{
  return "Delayed Group " + std::to_string(groups_[bin-1]);
}

//==============================================================================
// Fortran interoperability
//==============================================================================

extern "C" int delayedgroup_filter_groups(DelayedGroupFilter* filt, int i)
{return filt->groups_[i-1];}

} // namespace openmc
