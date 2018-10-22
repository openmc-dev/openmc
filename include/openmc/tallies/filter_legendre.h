#ifndef OPENMC_TALLIES_FILTER_LEGENDRE_H
#define OPENMC_TALLIES_FILTER_LEGENDRE_H

#include <string>

#include "openmc/cell.h"
#include "openmc/error.h"
#include "openmc/math_functions.h"
#include "openmc/tallies/filter.h"


namespace openmc {

//==============================================================================
//! Gives Legendre moments of the change in scattering angle
//==============================================================================

class LegendreFilter : public Filter
{
public:
  std::string type() const override {return "legendre";}

  ~LegendreFilter() = default;

  void
  from_xml(pugi::xml_node node) override
  {
    order_ = std::stoi(get_node_value(node, "order"));
    n_bins_ = order_ + 1;
  }

  void
  get_all_bins(Particle* p, int estimator, FilterMatch& match)
  const override
  {
    double wgt[n_bins_];
    calc_pn_c(order_, p->mu, wgt);
    for (int i = 0; i < n_bins_; i++) {
      match.bins_.push_back(i + 1);
      match.weights_.push_back(wgt[i]);
    }
  }

  void
  to_statepoint(hid_t filter_group) const override
  {
    Filter::to_statepoint(filter_group);
    write_dataset(filter_group, "order", order_);
  }

  std::string
  text_label(int bin) const override
  {
    return "Legendre expansion, P" + std::to_string(bin - 1);
  }

  int order_;
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_LEGENDRE_H
