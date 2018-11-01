#include "openmc/tallies/filter_sptl_legendre.h"

#include <utility>  // For pair

#include "openmc/capi.h"
#include "openmc/error.h"
#include "openmc/math_functions.h"
#include "openmc/xml_interface.h"

namespace openmc {

void
SpatialLegendreFilter::from_xml(pugi::xml_node node)
{
  order_ = std::stoi(get_node_value(node, "order"));

  auto axis = get_node_value(node, "axis");
  if (axis == "x") {
    axis_ = LegendreAxis::x;
  } else if (axis == "y") {
    axis_ = LegendreAxis::y;
  } else if (axis == "z") {
    axis_ = LegendreAxis::z;
  } else {
    fatal_error("Unrecognized axis on SpatialLegendreFilter");
  }

  min_ = std::stod(get_node_value(node, "min"));
  max_ = std::stod(get_node_value(node, "max"));

  n_bins_ = order_ + 1;
}

void
SpatialLegendreFilter::get_all_bins(const Particle* p, int estimator,
                                    FilterMatch& match) const
{
  // Get the coordinate along the axis of interest.
  double x;
  if (axis_ == LegendreAxis::x) {
    x = p->coord[0].xyz[0];
  } else if (axis_ == LegendreAxis::y) {
    x = p->coord[0].xyz[1];
  } else {
    x = p->coord[0].xyz[2];
  }

  if (x >= min_ && x <= max_) {
    // Compute the normalized coordinate value.
    double x_norm = 2.0*(x - min_) / (max_ - min_) - 1.0;

    // Compute and return the Legendre weights.
    double wgt[order_ + 1];
    calc_pn_c(order_, x_norm, wgt);
    for (int i = 0; i < order_ + 1; i++) {
      //TODO: off-by-one
      match.bins_.push_back(i + 1);
      match.weights_.push_back(wgt[i]);
    }
  }
}

void
SpatialLegendreFilter::to_statepoint(hid_t filter_group) const
{
  Filter::to_statepoint(filter_group);
  write_dataset(filter_group, "order", order_);
  if (axis_ == LegendreAxis::x) {
    write_dataset(filter_group, "axis", "x");
  } else if (axis_ == LegendreAxis::y) {
    write_dataset(filter_group, "axis", "y");
  } else {
    write_dataset(filter_group, "axis", "z");
  }
  write_dataset(filter_group, "min", min_);
  write_dataset(filter_group, "max", max_);
}

std::string
SpatialLegendreFilter::text_label(int bin) const
{
  std::stringstream out;
  out << "Legendre expansion, ";
  if (axis_ == LegendreAxis::x) {
    out << "x";
  } else if (axis_ == LegendreAxis::y) {
    out << "y";
  } else {
    out << "z";
  }
  //TODO: off-by-one
  out << " axis, P" << std::to_string(bin - 1);
  return out.str();
}

//==============================================================================
// C-API functions
//==============================================================================

std::pair<int, SpatialLegendreFilter*>
check_sptl_legendre_filter(int32_t index)
{
  // Make sure this is a valid index to an allocated filter.
  int err = verify_filter(index);
  if (err) {
    return {err, nullptr};
  }

  // Get a pointer to the filter and downcast.
  auto* filt_base = filter_from_f(index);
  auto* filt = dynamic_cast<SpatialLegendreFilter*>(filt_base);

  // Check the filter type.
  if (!filt) {
    set_errmsg("Not a spatial Legendre filter.");
    err = OPENMC_E_INVALID_TYPE;
  }
  return {err, filt};
}

extern "C" int
openmc_spatial_legendre_filter_get_order(int32_t index, int* order)
{
  // Check the filter.
  auto check_result = check_sptl_legendre_filter(index);
  auto err = check_result.first;
  auto filt = check_result.second;
  if (err) return err;

  // Output the order.
  *order = filt->order_;
  return 0;
}

extern "C" int
openmc_spatial_legendre_filter_get_params(int32_t index, int* axis,
                                          double* min, double* max)
{
  // Check the filter.
  auto check_result = check_sptl_legendre_filter(index);
  auto err = check_result.first;
  auto filt = check_result.second;
  if (err) return err;

  // Output the params.
  *axis = static_cast<int>(filt->axis_);
  *min = filt->min_;
  *max = filt->max_;
  return 0;
}

extern "C" int
openmc_spatial_legendre_filter_set_order(int32_t index, int order)
{
  // Check the filter.
  auto check_result = check_sptl_legendre_filter(index);
  auto err = check_result.first;
  auto filt = check_result.second;
  if (err) return err;

  // Update the filter.
  filt->order_ = order;
  filt->n_bins_ = order + 1;
  filter_update_n_bins(index);
  return 0;
}

extern "C" int
openmc_spatial_legendre_filter_set_params(int32_t index, const int* axis,
  const double* min, const double* max)
{
  // Check the filter.
  auto check_result = check_sptl_legendre_filter(index);
  auto err = check_result.first;
  auto filt = check_result.second;
  if (err) return err;

  // Update the filter.
  if (axis) filt->axis_ = static_cast<LegendreAxis>(*axis);
  if (min) filt->min_ = *min;
  if (max) filt->max_ = *max;
  return 0;
}

} // namespace openmc
