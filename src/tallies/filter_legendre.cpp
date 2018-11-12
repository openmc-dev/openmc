#include "openmc/tallies/filter_legendre.h"

#include "openmc/capi.h"
#include "openmc/error.h"
#include "openmc/math_functions.h"
#include "openmc/xml_interface.h"

namespace openmc {

void
LegendreFilter::from_xml(pugi::xml_node node)
{
  order_ = std::stoi(get_node_value(node, "order"));
  n_bins_ = order_ + 1;
}

void
LegendreFilter::get_all_bins(const Particle* p, int estimator,
                             FilterMatch& match) const
{
  double wgt[n_bins_];
  calc_pn_c(order_, p->mu, wgt);
  for (int i = 0; i < n_bins_; i++) {
    //TODO: off-by-one
    match.bins_.push_back(i + 1);
    match.weights_.push_back(wgt[i]);
  }
}

void
LegendreFilter::to_statepoint(hid_t filter_group) const
{
  Filter::to_statepoint(filter_group);
  write_dataset(filter_group, "order", order_);
}

std::string
LegendreFilter::text_label(int bin) const
{
  //TODO: off-by-one
  return "Legendre expansion, P" + std::to_string(bin - 1);
}

//==============================================================================
// C-API functions
//==============================================================================

extern "C" int
openmc_legendre_filter_get_order(int32_t index, int* order)
{
  // Make sure this is a valid index to an allocated filter.
  int err = verify_filter(index);
  if (err) return err;

  // Get a pointer to the filter and downcast.
  auto* filt_base = filter_from_f(index);
  auto* filt = dynamic_cast<LegendreFilter*>(filt_base);

  // Check the filter type.
  if (!filt) {
    set_errmsg("Not a legendre filter.");
    return OPENMC_E_INVALID_TYPE;
  }

  // Output the order.
  *order = filt->order_;
  return 0;
}

extern "C" int
openmc_legendre_filter_set_order(int32_t index, int order)
{
  // Make sure this is a valid index to an allocated filter.
  int err = verify_filter(index);
  if (err) return err;

  // Get a pointer to the filter and downcast.
  auto* filt_base = filter_from_f(index);
  auto* filt = dynamic_cast<LegendreFilter*>(filt_base);

  // Check the filter type.
  if (!filt) {
    set_errmsg("Not a legendre filter.");
    return OPENMC_E_INVALID_TYPE;
  }

  // Update the filter.
  filt->order_ = order;
  filt->n_bins_ = order + 1;
  filter_update_n_bins(index);
  return 0;
}

} // namespace openmc
