#ifndef OPENMC_TALLY_FILTER_SPTL_LEGENDRE_H
#define OPENMC_TALLY_FILTER_SPTL_LEGENDRE_H

#include <string>

#include "openmc/error.h"
#include "openmc/math_functions.h"
#include "openmc/tallies/tally_filter.h"


namespace openmc {

//TODO: those integer values are not needed when Fortran interop is removed
enum class LegendreAxis {
  x = 1, y = 2, z = 3
};

class SpatialLegendreFilter : public TallyFilter
{
public:
  virtual std::string type() const override {return "spatiallegendre";}

  virtual ~SpatialLegendreFilter() override = default;

  virtual void
  from_xml(pugi::xml_node node) override
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

  virtual void
  get_all_bins(Particle* p, int estimator, TallyFilterMatch& match)
  const override
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
        match.bins.push_back(i + 1);
        match.weights.push_back(wgt[i]);
      }
    }
  }

  virtual void
  to_statepoint(hid_t filter_group) const override
  {
    TallyFilter::to_statepoint(filter_group);
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

  virtual std::string
  text_label(int bin) const override
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

  int order_;
  LegendreAxis axis_;
  double min_, max_;
};

} // namespace openmc
#endif // OPENMC_TALLY_FILTER_SPTL_LEGENDRE_H
