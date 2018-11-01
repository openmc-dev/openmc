#include "openmc/tallies/filter_zernike.h"

#include <cmath>
#include <sstream>
#include <utility>  // For pair

#include "openmc/capi.h"
#include "openmc/error.h"
#include "openmc/math_functions.h"
#include "openmc/xml_interface.h"

namespace openmc {

//==============================================================================
// ZernikeFilter implementation
//==============================================================================

void
ZernikeFilter::from_xml(pugi::xml_node node)
{
  set_order(std::stoi(get_node_value(node, "order")));
  x_ = std::stod(get_node_value(node, "x"));
  y_ = std::stod(get_node_value(node, "y"));
  r_ = std::stod(get_node_value(node, "r"));
}

void
ZernikeFilter::get_all_bins(const Particle* p, int estimator,
                            FilterMatch& match) const
{
  // Determine the normalized (r,theta) coordinates.
  double x = p->coord[0].xyz[0] - x_;
  double y = p->coord[0].xyz[1] - y_;
  double r = std::sqrt(x*x + y*y) / r_;
  double theta = std::atan2(y, x);

  if (r <= 1.0) {
    // Compute and return the Zernike weights.
    double zn[n_bins_];
    calc_zn(order_, r, theta, zn);
    for (int i = 0; i < n_bins_; i++) {
      //TODO: off-by-one
      match.bins_.push_back(i+1);
      match.weights_.push_back(zn[i]);
    }
  }
}

void
ZernikeFilter::to_statepoint(hid_t filter_group) const
{
  Filter::to_statepoint(filter_group);
  write_dataset(filter_group, "order", order_);
  write_dataset(filter_group, "x", x_);
  write_dataset(filter_group, "y", y_);
  write_dataset(filter_group, "r", r_);
}

std::string
ZernikeFilter::text_label(int bin) const
{
  std::stringstream out;
  for (int n = 0; n < order_+1; n++) {
    int last = (n + 1) * (n + 2) / 2;
    //TODO: off-by-one
    if (bin <= last) {
      int first = last - n;
      int m = -n + (bin - first) * 2;
      out << "Zernike expansion, Z" << n << "," << m;
      return out.str();
    }
  }
}

void
ZernikeFilter::set_order(int order)
{
  order_ = order;
  n_bins_ = ((order+1) * (order+2)) / 2;
}

//==============================================================================
// ZernikeRadialFilter implementation
//==============================================================================

void
ZernikeRadialFilter::get_all_bins(const Particle* p, int estimator,
                                  FilterMatch& match) const
{
  // Determine the normalized radius coordinate.
  double x = p->coord[0].xyz[0] - x_;
  double y = p->coord[0].xyz[1] - y_;
  double r = std::sqrt(x*x + y*y) / r_;

  if (r <= 1.0) {
    // Compute and return the Zernike weights.
    double zn[n_bins_];
    calc_zn_rad(order_, r, zn);
    for (int i = 0; i < n_bins_; i++) {
      //TODO: off-by-one
      match.bins_.push_back(i+1);
      match.weights_.push_back(zn[i]);
    }
  }
}

std::string
ZernikeRadialFilter::text_label(int bin) const
{
  //TODO: off-by-one
  return "Zernike expansion, Z" + std::to_string(2*(bin-1)) + ",0";
}

void
ZernikeRadialFilter::set_order(int order)
{
  order_ = order;
  n_bins_ = order / 2 + 1;
}

//==============================================================================
// C-API functions
//==============================================================================

std::pair<int, ZernikeFilter*>
check_zernike_filter(int32_t index)
{
  // Make sure this is a valid index to an allocated filter.
  int err = verify_filter(index);
  if (err) {
    return {err, nullptr};
  }

  // Get a pointer to the filter and downcast.
  auto* filt_base = filter_from_f(index);
  auto* filt = dynamic_cast<ZernikeFilter*>(filt_base);

  // Check the filter type.
  if (!filt) {
    set_errmsg("Not a Zernike filter.");
    err = OPENMC_E_INVALID_TYPE;
  }
  return {err, filt};
}

extern "C" int
openmc_zernike_filter_get_order(int32_t index, int* order)
{
  // Check the filter.
  auto check_result = check_zernike_filter(index);
  auto err = check_result.first;
  auto filt = check_result.second;
  if (err) return err;

  // Output the order.
  *order = filt->order();
  return 0;
}

extern "C" int
openmc_zernike_filter_get_params(int32_t index, double* x, double* y,
                                 double* r)
{
  // Check the filter.
  auto check_result = check_zernike_filter(index);
  auto err = check_result.first;
  auto filt = check_result.second;
  if (err) return err;

  // Output the params.
  *x = filt->x_;
  *y = filt->y_;
  *r = filt->r_;
  return 0;
}

extern "C" int
openmc_zernike_filter_set_order(int32_t index, int order)
{
  // Check the filter.
  auto check_result = check_zernike_filter(index);
  auto err = check_result.first;
  auto filt = check_result.second;
  if (err) return err;

  // Update the filter.
  filt->set_order(order);
  filter_update_n_bins(index);
  return 0;
}

extern "C" int
openmc_zernike_filter_set_params(int32_t index, const double* x,
                                 const double* y, const double* r)
{
  // Check the filter.
  auto check_result = check_zernike_filter(index);
  auto err = check_result.first;
  auto filt = check_result.second;
  if (err) return err;

  // Update the filter.
  if (x) filt->x_ = *x;
  if (y) filt->y_ = *y;
  if (r) filt->r_ = *r;
  return 0;
}

} // namespace openmc
