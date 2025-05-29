#include "openmc/tallies/filter_weight.h"

#include <algorithm> // for is_sorted
#include <stdexcept> // for runtime_error

#include <fmt/core.h>

#include "openmc/search.h"
#include "openmc/xml_interface.h"

namespace openmc {

//==============================================================================
// WeightFilter implementation
//==============================================================================

void WeightFilter::from_xml(pugi::xml_node node)
{
  auto bins = get_node_array<double>(node, "bins");
  this->set_bins(bins);
}

void WeightFilter::set_bins(span<const double> bins)
{
  if (!std::is_sorted(bins.begin(), bins.end())) {
    throw std::runtime_error {"Weight bins must be monotonically increasing."};
  }

  // Clear existing bins
  bins_.clear();
  bins_.reserve(bins.size());

  // Copy bins
  bins_.insert(bins_.end(), bins.begin(), bins.end());
  n_bins_ = bins_.size() - 1;
}

void WeightFilter::get_all_bins(
  const Particle& p, TallyEstimator estimator, FilterMatch& match) const
{
  // Get particle weight
  double wgt = p.wgt_last();

  // Bin the weight
  if (wgt >= bins_.front() && wgt <= bins_.back()) {
    auto bin = lower_bound_index(bins_.begin(), bins_.end(), wgt);
    match.bins_.push_back(bin);
    match.weights_.push_back(1.0);
  }
}

void WeightFilter::to_statepoint(hid_t filter_group) const
{
  Filter::to_statepoint(filter_group);
  write_dataset(filter_group, "bins", bins_);
}

std::string WeightFilter::text_label(int bin) const
{
  return fmt::format("Weight [{}, {}]", bins_[bin], bins_[bin + 1]);
}

} // namespace openmc
