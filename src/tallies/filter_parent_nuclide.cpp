#include "openmc/tallies/filter_parent_nuclide.h"

#include <fmt/core.h>

#include "openmc/capi.h"
#include "openmc/chain.h"
#include "openmc/search.h"
#include "openmc/settings.h"
#include "openmc/xml_interface.h"

namespace openmc {

//==============================================================================
// ParentNuclideFilter implementation
//==============================================================================

void ParentNuclideFilter::from_xml(pugi::xml_node node)
{
  nuclides_ = get_node_array<std::string>(node, "bins");

  // Convert nuclides to indices in data::chain_nuclides
  std::vector<int> bins;
  for (const auto& nuclide : nuclides_) {
    auto it = data::chain_nuclide_map.find(nuclide);
    if (it != data::chain_nuclide_map.end()) {
      bins.push_back(it->second);
    } else {
      // The default value of parent_nuclide is -1, so to prevent a score to
      // this bin assign the value -2.
      bins.push_back(-2);
    }
  }
  this->set_bins(bins);
}

void ParentNuclideFilter::set_bins(gsl::span<const int> bins)
{
  // Clear existing bins
  bins_.clear();
  bins_.reserve(bins.size());
  map_.clear();

  // Set bins based on chain nuclide indexing
  for (gsl::index i = 0; i < bins.size(); ++i) {
    bins_.push_back(bins[i]);
    map_[bins[i]] = i;
  }

  n_bins_ = bins_.size();
}

void ParentNuclideFilter::get_all_bins(
  const Particle& p, TallyEstimator estimator, FilterMatch& match) const
{
  // Get the particle's parent nuclide
  int parent_nuclide = p.parent_nuclide();

  // Find bin matching parent nuclide
  auto search = map_.find(parent_nuclide);
  if (search != map_.end()) {
    match.bins_.push_back(search->second);
    match.weights_.push_back(1.0);
  }
}

void ParentNuclideFilter::to_statepoint(hid_t filter_group) const
{
  Filter::to_statepoint(filter_group);
  write_dataset(filter_group, "bins", nuclides_);
}

std::string ParentNuclideFilter::text_label(int bin) const
{
  return fmt::format("Parent Nuclide {}", nuclides_[bin]);
}

} // namespace openmc
