#include "openmc/tallies/filter_time.h"

#include <algorithm>  // for min, max, copy, adjacent_find
#include <functional> // for greater_equal
#include <iterator>   // for back_inserter

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

  // Ensure time bins are sorted and don't have duplicates
  if (std::adjacent_find(bins.cbegin(), bins.cend(), std::greater_equal<>()) !=
      bins.end()) {
    throw std::runtime_error {"Time bins must be monotonically increasing."};
  }

  // Copy bins
  std::copy(bins.cbegin(), bins.cend(), std::back_inserter(bins_));
  n_bins_ = bins_.size() - 1;
}

void TimeFilter::get_all_bins(
  const Particle& p, TallyEstimator estimator, FilterMatch& match) const
{
  // Get the start/end time of the particle for this track
  const auto t_start = p.time_last();
  const auto t_end = p.time();

  // If time interval is entirely out of time bin range, exit
  if (t_end < bins_.front() || t_start >= bins_.back())
    return;

  if (estimator == TallyEstimator::TRACKLENGTH) {
    // -------------------------------------------------------------------------
    // For tracklength estimator, we have to check the start/end time of
    // the current track and find where it overlaps with time bins and score
    // accordingly

    // Determine first bin containing a portion of time interval
    auto i_bin = lower_bound_index(bins_.begin(), bins_.end(), t_start);

    // If time interval is zero, add a match corresponding to the starting time
    if (t_end == t_start) {
      match.set(i_bin);
      return;
    }

    /* TODO: this is incorrect when used in conjunction with other
     * vector-matching that have spatial information. The problem is
     * that the spatial matches are correlated with time, so we
     * are in fact dealing with correlated matches across multiple
     * filters: something that OpenMC cannot and should not handle.
     *
     * With, for example, an energy functionaly expansion filter,
     * this is correct, because the vector match to the energy is
     * not changing over time (assuming we're dealing with neutral
     * particles)
     */

    // Find matching bins. There are multiple matching bins here,
    // so we use the vector interface on FilterMatch.
    double dt_total = t_end - t_start;
    auto& match_vector = match.vector_pairs();
    for (; i_bin < bins_.size() - 1; ++i_bin) {
      double t_left = std::max(t_start, bins_[i_bin]);
      double t_right = std::min(t_end, bins_[i_bin + 1]);

      // Add match with weight equal to the fraction of the time interval within
      // the current time bin
      const double fraction = (t_right - t_left) / dt_total;
      match_vector.push_back({i_bin, fraction});

      if (t_end < bins_[i_bin + 1])
        break;
    }
  } else if (t_end < bins_.back()) {
    // -------------------------------------------------------------------------
    // For collision estimator or surface tallies, find a match based on the
    // exact time of the particle

    const auto i_bin = lower_bound_index(bins_.begin(), bins_.end(), t_end);
    match.set(i_bin);
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
