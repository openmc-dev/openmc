#include "openmc/tallies/filter.h"

#include "openmc/capi.h"
#include "openmc/error.h"
#include "openmc/math_functions.h"
#include "openmc/xml_interface.h"

namespace openmc {

void
Filter::LegendreFilter_from_xml(pugi::xml_node node)
{
  this->set_order(std::stoi(get_node_value(node, "order")));
}

void
Filter::LegendreFilter_set_order(int order)
{
  if (order < 0) {
    throw std::invalid_argument{"Legendre order must be non-negative."};
  }
  order_ = order;
  n_bins_ = order_ + 1;
}

void
Filter::LegendreFilter_get_all_bins(const Particle& p, TallyEstimator estimator,
                             FilterMatch& match) const
{
  double wgt[FILTERMATCH_BINS_WEIGHTS_SIZE];
  calc_pn_c(order_, p.mu_, wgt);
  for (int i = 0; i < n_bins_; i++) {
    //match.bins_.push_back(i);
    //match.weights_.push_back(wgt[i]);
    match.push_back(i, wgt[i]);
  }
}

void
Filter::LegendreFilter_to_statepoint(hid_t filter_group) const
{
  write_dataset(filter_group, "order", order_);
}

std::string
Filter::LegendreFilter_text_label(int bin) const
{
  return "Legendre expansion, P" + std::to_string(bin);
}

//==============================================================================
// C-API functions
//==============================================================================

/*

extern "C" int
openmc_legendre_filter_get_order(int32_t index, int* order)
{
  // Make sure this is a valid index to an allocated filter.
  if (int err = verify_filter(index)) return err;

  // Get a pointer to the filter and downcast.
  const auto& filt_base = model::tally_filters[index].get();
  auto* filt = dynamic_cast<LegendreFilter*>(filt_base);

  // Check the filter type.
  if (!filt) {
    set_errmsg("Not a legendre filter.");
    return OPENMC_E_INVALID_TYPE;
  }

  // Output the order.
  *order = filt->order();
  return 0;
}

extern "C" int
openmc_legendre_filter_set_order(int32_t index, int order)
{
  // Make sure this is a valid index to an allocated filter.
  if (int err = verify_filter(index)) return err;

  // Get a pointer to the filter and downcast.
  const auto& filt_base = model::tally_filters[index].get();
  auto* filt = dynamic_cast<LegendreFilter*>(filt_base);

  // Check the filter type.
  if (!filt) {
    set_errmsg("Not a legendre filter.");
    return OPENMC_E_INVALID_TYPE;
  }

  // Update the filter.
  filt->set_order(order);
  return 0;
}
*/

} // namespace openmc
