#ifndef OPENMC_TALLY_FILTER_ZERNIKE_H
#define OPENMC_TALLY_FILTER_ZERNIKE_H

#include <cmath>
#include <sstream>
#include <string>

#include "openmc/error.h"
#include "openmc/math_functions.h"
#include "openmc/tallies/tally_filter.h"


namespace openmc {

class ZernikeFilter : public TallyFilter
{
public:
  virtual std::string type() const override {return "zernike";}

  virtual ~ZernikeFilter() override = default;

  virtual void
  from_xml(pugi::xml_node node) override
  {
    order_ = std::stoi(get_node_value(node, "order"));
    x_ = std::stod(get_node_value(node, "x"));
    y_ = std::stod(get_node_value(node, "y"));
    r_ = std::stod(get_node_value(node, "r"));
    calc_n_bins();
  }

  virtual void
  get_all_bins(Particle* p, int estimator, TallyFilterMatch& match)
  const override
  {
    // Determine the normalized (r,theta) coordinates.
    double x = p->coord[0].xyz[0] - x_;
    double y = p->coord[0].xyz[1] - y_;
    double r = std::sqrt(x*x + y*y) / r_;
    double theta = std::atan2(y, x);

    if (r <= 1.0) {
      // Compute and return the Zernike weights.
      double zn[n_bins_];
      calc_zn_c(order_, r, theta, zn);
      for (int i = 0; i < n_bins_; i++) {
        match.bins.push_back(i+1);
        match.weights.push_back(zn[i]);
      }
    }
  }

  virtual void
  to_statepoint(hid_t filter_group) const override
  {
    TallyFilter::to_statepoint(filter_group);
    write_dataset(filter_group, "order", order_);
    write_dataset(filter_group, "x", x_);
    write_dataset(filter_group, "y", y_);
    write_dataset(filter_group, "r", r_);
  }

  virtual std::string
  text_label(int bin) const override
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

  virtual void calc_n_bins() {n_bins_ = ((order_+1) * (order_+2)) / 2;}

  int order_;
  double x_, y_, r_;
};

class ZernikeRadialFilter : public ZernikeFilter
{
public:
  virtual std::string type() const override {return "zernikeradial";}

  virtual ~ZernikeRadialFilter() override = default;

  virtual void
  get_all_bins(Particle* p, int estimator, TallyFilterMatch& match)
  const override
  {
    // Determine the normalized radius coordinate.
    double x = p->coord[0].xyz[0] - x_;
    double y = p->coord[0].xyz[1] - y_;
    double r = std::sqrt(x*x + y*y) / r_;

    if (r <= 1.0) {
      // Compute and return the Zernike weights.
      double zn[n_bins_];
      calc_zn_rad_c(order_, r, zn);
      for (int i = 0; i < n_bins_; i++) {
        match.bins.push_back(i+1);
        match.weights.push_back(zn[i]);
      }
    }
  }

  virtual std::string
  text_label(int bin) const override
  {
    return "Zernike expansion, Z" + std::to_string(2*(bin-1)) + ",0";
  }

  virtual void calc_n_bins() override {n_bins_ = order_ / 2 + 1;}
};

} // namespace openmc
#endif // OPENMC_TALLY_FILTER_ZERNIKE_H
