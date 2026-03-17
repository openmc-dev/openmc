#include "openmc/tallies/filter_sptl_fourier.h"

#include <utility> // For pair

#include <fmt/core.h>

#include "openmc/capi.h"
#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/xml_interface.h"

namespace openmc {

void SpatialFourierFilter::from_xml(pugi::xml_node node)
{
  this->set_order(std::stoi(get_node_value(node, "order")));

  auto axis = get_node_value(node, "axis");
  switch (axis[0]) {
  case 'x':
    this->set_axis(FourierAxis::x);
    break;
  case 'y':
    this->set_axis(FourierAxis::y);
    break;
  case 'z':
    this->set_axis(FourierAxis::z);
    break;
  default:
    throw std::runtime_error {
      "Axis for SpatialFourierFilter must be 'x', 'y', or 'z'"};
  }

  double min = std::stod(get_node_value(node, "min"));
  double max = std::stod(get_node_value(node, "max"));
  this->set_minmax(min, max);
}

void SpatialFourierFilter::set_order(int order)
{
  if (order < 0) {
    throw std::invalid_argument {"Fourier order must be non-negative."};
  }
  order_ = order;
  n_bins_ = 2 * order_ + 1;
}

void SpatialFourierFilter::set_axis(FourierAxis axis)
{
  axis_ = axis;
}

void SpatialFourierFilter::set_minmax(double min, double max)
{
  if (max <= min) {
    throw std::invalid_argument {
      "Maximum value must be greater than minimum value"};
  }
  min_ = min;
  max_ = max;
}

void SpatialFourierFilter::get_all_bins(
  const Particle& p, TallyEstimator estimator, FilterMatch& match) const
{
  // Get the coordinate along the axis of interest.
  double x;
  if (axis_ == FourierAxis::x) {
    x = p.r().x;
  } else if (axis_ == FourierAxis::y) {
    x = p.r().y;
  } else {
    x = p.r().z;
  }

  if (x >= min_ && x <= max_) {
    // Compute the normalized coordinate value on [0, 1]
    double x_norm = (x - min_) / (max_ - min_);

    // Compute and return the Fourier weights.
    vector<double> wgt(n_bins_);
    wgt[0] = 1.0;  // a_0: constant term
    for (int n = 1; n <= order_; ++n) {
      double arg = 2.0 * PI * n * x_norm;
      wgt[2*n - 1] = std::cos(arg);
      wgt[2*n] = std::sin(arg);
    }
    for (int i = 0; i < n_bins_; ++i) {
      match.bins_.push_back(i);
      match.weights_.push_back(wgt[i]);
    }
  }
}

void SpatialFourierFilter::to_statepoint(hid_t filter_group) const
{
  Filter::to_statepoint(filter_group);
  write_dataset(filter_group, "order", order_);
  if (axis_ == FourierAxis::x) {
    write_dataset(filter_group, "axis", "x");
  } else if (axis_ == FourierAxis::y) {
    write_dataset(filter_group, "axis", "y");
  } else {
    write_dataset(filter_group, "axis", "z");
  }
  write_dataset(filter_group, "min", min_);
  write_dataset(filter_group, "max", max_);
}

std::string SpatialFourierFilter::text_label(int bin) const
{
  std::string axis_str;
  std::string func_str;
  if (axis_ == FourierAxis::x) {
    axis_str = "x";
  } else if (axis_ == FourierAxis::y) {
    axis_str = "y";
  } else {
    axis_str = "z";
  }

  if (bin == 0) {
    func_str = "a0 (constant)";
  } else if (bin % 2 == 1) {
    int n = (bin + 1) / 2;
    func_str = fmt::format("a{} (cos)", n);
  } else {
    int n = bin / 2;
    func_str = fmt::format("b{} (sin)", n);
  }
  return fmt::format("Fourier expansion, {} axis, {}", axis_str, func_str);
}

//==============================================================================
// C-API functions
//==============================================================================

std::pair<int, SpatialFourierFilter*> check_sptl_fourier_filter(int32_t index)
{
  // Make sure this is a valid index to an allocated filter.
  int err = verify_filter(index);
  if (err) {
    return {err, nullptr};
  }

  // Get a pointer to the filter and downcast.
  const auto& filt_base = model::tally_filters[index].get();
  auto* filt = dynamic_cast<SpatialFourierFilter*>(filt_base);

  // Check the filter type.
  if (!filt) {
    set_errmsg("Not a spatial Fourier filter.");
    err = OPENMC_E_INVALID_TYPE;
  }
  return {err, filt};
}

extern "C" int openmc_spatial_fourier_filter_get_order(
  int32_t index, int* order)
{
  // Check the filter.
  auto check_result = check_sptl_fourier_filter(index);
  auto err = check_result.first;
  auto filt = check_result.second;
  if (err)
    return err;

  // Output the order.
  *order = filt->order();
  return 0;
}

extern "C" int openmc_spatial_fourier_filter_get_params(
  int32_t index, int* axis, double* min, double* max)
{
  // Check the filter.
  auto check_result = check_sptl_fourier_filter(index);
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

extern "C" int openmc_spatial_fourier_filter_set_order(
  int32_t index, int order)
{
  // Check the filter.
  auto check_result = check_sptl_fourier_filter(index);
  auto err = check_result.first;
  auto filt = check_result.second;
  if (err)
    return err;

  // Update the filter.
  filt->set_order(order);
  return 0;
}

extern "C" int openmc_spatial_fourier_filter_set_params(
  int32_t index, const int* axis, const double* min, const double* max)
{
  // Check the filter.
  auto check_result = check_sptl_fourier_filter(index);
  auto err = check_result.first;
  auto filt = check_result.second;
  if (err)
    return err;

  // Update the filter.
  if (axis)
    filt->set_axis(static_cast<FourierAxis>(*axis));
  if (min && max)
    filt->set_minmax(*min, *max);
  return 0;
}

} // namespace openmc
