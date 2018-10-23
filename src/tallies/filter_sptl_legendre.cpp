#include "openmc/tallies/filter_sptl_legendre.h"

#include "openmc/capi.h"
#include "openmc/error.h"
#include "openmc/math_functions.h"

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
SpatialLegendreFilter::get_all_bins(Particle* p, int estimator, FilterMatch& match)
const
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
  out << " axis, P" << std::to_string(bin - 1);
  return out.str();
}

//==============================================================================
// C-API functions
//==============================================================================

extern "C" int
openmc_spatial_legendre_filter_get_order(int32_t index, int* order)
{
  int err = verify_filter(index);
  if (err) return err;

  auto filt = filter_from_f(index);
  if (filt->type() != "spatiallegendre") {
    set_errmsg("Not a spatial Legendre filter.");
    return OPENMC_E_INVALID_TYPE;
  }

  auto l_filt = static_cast<SpatialLegendreFilter*>(filt);
  *order = l_filt->order_;
  return 0;
}

extern "C" int
openmc_spatial_legendre_filter_get_params(int32_t index, int* axis,
                                          double* min, double* max)
{
  int err = verify_filter(index);
  if (err) return err;

  auto filt = filter_from_f(index);
  if (filt->type() != "spatiallegendre") {
    set_errmsg("Not a spatial Legendre filter.");
    return OPENMC_E_INVALID_TYPE;
  }

  auto l_filt = static_cast<SpatialLegendreFilter*>(filt);
  *axis = static_cast<int>(l_filt->axis_);
  *min = l_filt->min_;
  *max = l_filt->max_;
  return 0;
}

extern "C" int
openmc_spatial_legendre_filter_set_order(int32_t index, int order)
{
  int err = verify_filter(index);
  if (err) return err;

  auto filt = filter_from_f(index);
  if (filt->type() != "spatiallegendre") {
    set_errmsg("Not a spatial Legendre filter.");
    return OPENMC_E_INVALID_TYPE;
  }

  auto l_filt = static_cast<SpatialLegendreFilter*>(filt);
  l_filt->order_ = order;
  l_filt->n_bins_ = order + 1;
  filter_update_n_bins(index);
  return 0;
}

extern "C" int
openmc_spatial_legendre_filter_set_params(int32_t index, const int* axis,
  const double* min, const double* max)
{
  int err = verify_filter(index);
  if (err) return err;

  auto filt = filter_from_f(index);
  if (filt->type() != "spatiallegendre") {
    set_errmsg("Not a spatial Legendre filter.");
    return OPENMC_E_INVALID_TYPE;
  }

  auto l_filt = static_cast<SpatialLegendreFilter*>(filt);
  if (axis) l_filt->axis_ = static_cast<LegendreAxis>(*axis);
  if (min) l_filt->min_ = *min;
  if (max) l_filt->max_ = *max;
  return 0;
}

} // namespace openmc
