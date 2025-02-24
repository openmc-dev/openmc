#include "openmc/tallies/filter_weight.h"

#include <fmt/core.h>

#include "openmc/capi.h"
#include "openmc/search.h"
#include "openmc/constants.h"
#include "openmc/settings.h"
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

void WeightFilter::set_bins(gsl::span<const double> bins)
{
   // Clear existing bins
   bins_.clear();
   bins_.reserve(bins.size());
 
   // Copy bins
   for (gsl::index i = 0; i < bins.size(); ++i) {
     bins_.push_back(bins[i]);
   }
 
   n_bins_ = bins_.size();
}

void WeightFilter::get_all_bins(const Particle& p, TallyEstimator estimator, FilterMatch& match) const
{
  // Get particle weight
  double wgt = p.wgt_last(); // or wgt_last()?

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
  return fmt::format("Weight bin [{}, {}]", bins_[bin], bins_[bin+1]);
}

} // namespace openmc
