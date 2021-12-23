#include "openmc/tallies/filter_time.h"

#include <algorithm> // for min, max

#include <fmt/core.h>

#include "openmc/search.h"
#include "openmc/xml_interface.h"

namespace openmc {

//==============================================================================
// TimeFilter implementation
//==============================================================================

void TimeFilter::from_xml(pugi::xml_node node)
{
  auto bins = get_node_array<double>(node, "bins");
  this->set_bins(bins);
}

void TimeFilter::set_bins(gsl::span<const double> bins)
{
  // Clear existing bins
  bins_.clear();
  bins_.reserve(bins.size());

  // Copy bins, ensuring they are valid
  for (gsl::index i = 0; i < bins.size(); ++i) {
    if (i > 0 && bins[i] <= bins[i - 1]) {
      throw std::runtime_error {"Time bins must be monotonically increasing."};
    }
    bins_.push_back(bins[i]);
  }

  n_bins_ = bins_.size() - 1;
}

void TimeFilter::get_all_bins(
  const Particle& p, TallyEstimator estimator, FilterMatch& match) const
{
  // Get the start/end time of the particle for this track
  auto t_start = p.time_last();
  auto t_end = p.time();

  // If time interval is entirely out of time bin range, exit
  if (t_end < bins_.front() || t_start >= bins_.back())
    return;

  if (estimator == TallyEstimator::ANALOG) {
    // -------------------------------------------------------------------------
    // For surface tallies, find a match based on the exact time the particle
    // crosses the surface
    auto i_bin = lower_bound_index(bins_.begin(), bins_.end(), t_end);
    match.bins_.push_back(i_bin);
    match.weights_.push_back(1.0);

  } else {
    // -------------------------------------------------------------------------
    // For volume tallies, we have to check the start/end time of the current
    // track and find where it overlaps with time bins and score accordingly

    // Determine first bin containing a portion of time interval
    auto i_bin = lower_bound_index(bins_.begin(), bins_.end(), t_start);

    // Find matching bins
    double dt_total = t_end - t_start;
    for (; i_bin < bins_.size() - 1; ++i_bin) {
      double t_left = std::max(t_start, bins_[i_bin]);
      double t_right = std::min(t_end, bins_[i_bin + 1]);

      // Add match with weight equal to the fraction of the time interval within
      // the current time bin
      if (dt_total > 0.0) {
        double fraction = (t_right - t_left) / dt_total;
        match.bins_.push_back(i_bin);
        match.weights_.push_back(fraction);
      }

      if (t_end < bins_[i_bin + 1])
        break;
    }
  }
}

void TimeFilter::to_statepoint(hid_t filter_group) const
{
  Filter::to_statepoint(filter_group);
  write_dataset(filter_group, "bins", bins_);
}

std::string TimeFilter::text_label(int bin) const
{
  return fmt::format("Time [{}, {})", bins_[bin], bins_[bin + 1]);
}

} // namespace openmc
