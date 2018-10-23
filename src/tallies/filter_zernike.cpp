#include "openmc/tallies/filter_zernike.h"

#include <cmath>
#include <sstream>

#include "openmc/capi.h"
#include "openmc/error.h"
#include "openmc/math_functions.h"

namespace openmc {

//==============================================================================
// ZernikeFilter implementation
//==============================================================================

void
ZernikeFilter::from_xml(pugi::xml_node node)
{
  order_ = std::stoi(get_node_value(node, "order"));
  x_ = std::stod(get_node_value(node, "x"));
  y_ = std::stod(get_node_value(node, "y"));
  r_ = std::stod(get_node_value(node, "r"));
  calc_n_bins();
}

void
ZernikeFilter::get_all_bins(Particle* p, int estimator, FilterMatch& match)
const
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
    if (bin <= last) {
      int first = last - n;
      int m = -n + (bin - first) * 2;
      out << "Zernike expansion, Z" << n << "," << m;
      return out.str();
    }
  }
}

//==============================================================================
// ZernikeRadialFilter implementation
//==============================================================================

void
ZernikeRadialFilter::get_all_bins(Particle* p, int estimator,
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
      match.bins_.push_back(i+1);
      match.weights_.push_back(zn[i]);
    }
  }
}

std::string
ZernikeRadialFilter::text_label(int bin) const
{
  return "Zernike expansion, Z" + std::to_string(2*(bin-1)) + ",0";
}

//==============================================================================
// C-API functions
//==============================================================================

extern "C" int
openmc_zernike_filter_get_order(int32_t index, int* order)
{
  int err = verify_filter(index);
  if (err) return err;

  auto filt = filter_from_f(index);
  if (filt->type() != "zernike" && filt->type() != "zernikeradial") {
    set_errmsg("Not a Zernike filter.");
    return OPENMC_E_INVALID_TYPE;
  }

  auto z_filt = static_cast<ZernikeFilter*>(filt);
  *order = z_filt->order_;
  return 0;
}

extern "C" int
openmc_zernike_filter_get_params(int32_t index, double* x, double* y,
                                 double* r)
{
  int err = verify_filter(index);
  if (err) return err;

  auto filt = filter_from_f(index);
  if (filt->type() != "zernike" && filt->type() != "zernikeradial") {
    set_errmsg("Not a Zernike filter.");
    return OPENMC_E_INVALID_TYPE;
  }

  auto z_filt = static_cast<ZernikeFilter*>(filt);
  *x = z_filt->x_;
  *y = z_filt->y_;
  *r = z_filt->r_;
  return 0;
}

extern "C" int
openmc_zernike_filter_set_order(int32_t index, int order)
{
  int err = verify_filter(index);
  if (err) return err;

  auto filt = filter_from_f(index);
  if (filt->type() != "zernike" && filt->type() != "zernikeradial") {
    set_errmsg("Not a Zernike filter.");
    return OPENMC_E_INVALID_TYPE;
  }

  auto z_filt = static_cast<ZernikeFilter*>(filt);
  z_filt->order_ = order;
  z_filt->calc_n_bins();
  filter_update_n_bins(index);
  return 0;
}

extern "C" int
openmc_zernike_filter_set_params(int32_t index, const double* x,
                                 const double* y, const double* r)
{
  int err = verify_filter(index);
  if (err) return err;

  auto filt = filter_from_f(index);
  if (filt->type() != "zernike" && filt->type() != "zernikeradial") {
    set_errmsg("Not a Zernike filter.");
    return OPENMC_E_INVALID_TYPE;
  }

  auto z_filt = static_cast<ZernikeFilter*>(filt);
  if (x) z_filt->x_ = *x;
  if (y) z_filt->y_ = *y;
  if (r) z_filt->r_ = *r;
  return 0;
}

} // namespace openmc
