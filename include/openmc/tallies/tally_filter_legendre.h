#ifndef OPENMC_TALLY_FILTER_LEGENDRE_H
#define OPENMC_TALLY_FILTER_LEGENDRE_H

#include <string>

#include "openmc/cell.h"
#include "openmc/error.h"
#include "openmc/math_functions.h"
#include "openmc/tallies/tally_filter.h"


namespace openmc {

//==============================================================================
//! Gives Legendre moments of the change in scattering angle
//==============================================================================

class LegendreFilter : public TallyFilter
{
public:
  virtual std::string type() const override {return "legendre";}

  virtual ~LegendreFilter() override = default;

  virtual void
  from_xml(pugi::xml_node node) override
  {
    order_ = std::stoi(get_node_value(node, "order"));
    n_bins_ = order_ + 1;
  }

  virtual void
  get_all_bins(Particle* p, int estimator, TallyFilterMatch& match)
  const override
  {
    double wgt[n_bins_];
    calc_pn_c(order_, p->mu, wgt);
    for (int i = 0; i < n_bins_; i++) {
      match.bins_.push_back(i + 1);
      match.weights_.push_back(wgt[i]);
    }
  }

  virtual void
  to_statepoint(hid_t filter_group) const override
  {
    TallyFilter::to_statepoint(filter_group);
    write_dataset(filter_group, "order", order_);
  }

  virtual std::string
  text_label(int bin) const override
  {
    return "Legendre expansion, P" + std::to_string(bin - 1);
  }

  int order_;
};

} // namespace openmc
#endif // OPENMC_TALLY_FILTER_LEGENDRE_H
