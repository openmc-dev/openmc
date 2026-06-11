#include "openmc/tallies/filter_sptl_legendre.h"

#include <utility> // For pair

#include <fmt/core.h>

#include "openmc/capi.h"
#include "openmc/error.h"
#include "openmc/math_functions.h"

namespace openmc {

void SpatialLegendreFilter::set_order(int order)
{
  if (order < 0) {
    throw std::invalid_argument {"Legendre order must be non-negative."};
  }
  order_ = order;
  n_bins_ = order_ + 1;
}

void SpatialLegendreFilter::get_all_bins(
  const Particle& p, TallyEstimator estimator, FilterMatch& match) const
{
  // Get the coordinate along the axis of interest.
  double x = this->position(p);

  if (x >= min_ && x <= max_) {
    // Compute the normalized coordinate value.
    double x_norm = 2.0 * (x - min_) / (max_ - min_) - 1.0;

    // Compute and return the Legendre weights.
    vector<double> wgt(order_ + 1);
    calc_pn_c(order_, x_norm, wgt.data());
    for (int i = 0; i < order_ + 1; i++) {
      match.bins_.push_back(i);
      match.weights_.push_back(wgt[i]);
    }
  }
}

std::string SpatialLegendreFilter::text_label(int bin) const
{
  return fmt::format(
    "Legendre expansion, {} axis, P{}", this->axis_label(), bin);
}

//==============================================================================
// C-API functions
//==============================================================================

std::pair<int, SpatialLegendreFilter*> check_sptl_legendre_filter(int32_t index)
{
  // Make sure this is a valid index to an allocated filter.
  int err = verify_filter(index);
  if (err) {
    return {err, nullptr};
  }

  // Get a pointer to the filter and downcast.
  const auto& filt_base = model::tally_filters[index].get();
  auto* filt = dynamic_cast<SpatialLegendreFilter*>(filt_base);

  // Check the filter type.
  if (!filt) {
    set_errmsg("Not a spatial Legendre filter.");
    err = OPENMC_E_INVALID_TYPE;
  }
  return {err, filt};
}

extern "C" int openmc_spatial_legendre_filter_get_order(
  int32_t index, int* order)
{
  // Check the filter.
  auto check_result = check_sptl_legendre_filter(index);
  auto err = check_result.first;
  auto filt = check_result.second;
  if (err)
    return err;

  // Output the order.
  *order = filt->order();
  return 0;
}

extern "C" int openmc_spatial_legendre_filter_get_params(
  int32_t index, int* axis, double* min, double* max)
{
  // Check the filter.
  auto check_result = check_sptl_legendre_filter(index);
  auto err = check_result.first;
  auto filt = check_result.second;
  if (err)
    return err;

  // Output the params.
  *axis = static_cast<int>(filt->axis());
  *min = filt->min();
  *max = filt->max();
  return 0;
}

extern "C" int openmc_spatial_legendre_filter_set_order(
  int32_t index, int order)
{
  // Check the filter.
  auto check_result = check_sptl_legendre_filter(index);
  auto err = check_result.first;
  auto filt = check_result.second;
  if (err)
    return err;

  // Update the filter.
  filt->set_order(order);
  return 0;
}

extern "C" int openmc_spatial_legendre_filter_set_params(
  int32_t index, const int* axis, const double* min, const double* max)
{
  // Check the filter.
  auto check_result = check_sptl_legendre_filter(index);
  auto err = check_result.first;
  auto filt = check_result.second;
  if (err)
    return err;

  // Update the filter.
  if (axis)
    filt->set_axis(static_cast<SpatialExpansionFilter::Axis>(*axis));
  if (min && max)
    filt->set_minmax(*min, *max);
  return 0;
}

} // namespace openmc
