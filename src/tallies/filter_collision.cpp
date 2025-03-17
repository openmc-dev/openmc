#include "openmc/tallies/filter_collision.h"

#include <fmt/core.h>

#include "openmc/capi.h"
#include "openmc/search.h"
#include "openmc/settings.h"
#include "openmc/xml_interface.h"

namespace openmc {

//==============================================================================
// CollisionFilter implementation
//==============================================================================

void CollisionFilter::from_xml(pugi::xml_node node)
{
  auto bins = get_node_array<int>(node, "bins");
  this->set_bins(bins);
}

void CollisionFilter::set_bins(span<const int> bins)
{
  // Clear existing bins
  bins_.clear();
  bins_.reserve(bins.size());
  map_.clear();

  // Copy bins
  for (int64_t i = 0; i < bins.size(); ++i) {
    bins_.push_back(bins[i]);
    map_[bins[i]] = i;
  }

  n_bins_ = bins_.size();
}

void CollisionFilter::get_all_bins(
  const Particle& p, TallyEstimator estimator, FilterMatch& match) const
{
  // Get the number of collisions for the particle
  auto n = p.n_collision();

  // Bin the collision number. Must fit exactly the desired collision number.
  auto search = map_.find(n);
  if (search != map_.end()) {
    match.bins_.push_back(search->second);
    match.weights_.push_back(1.0);
  }
}

void CollisionFilter::to_statepoint(hid_t filter_group) const
{
  Filter::to_statepoint(filter_group);
  write_dataset(filter_group, "bins", bins_);
}

std::string CollisionFilter::text_label(int bin) const
{
  return fmt::format("Collision Number {}", bins_[bin]);
}

} // namespace openmc
