#include "openmc/tallies/filter_lifetime_moment.h"

#include <stdexcept> // for invalid_argument

#include "openmc/capi.h"
#include "openmc/error.h"
#include "openmc/xml_interface.h"

namespace openmc {

//==============================================================================
// LifetimeMomentFilter implementation
//==============================================================================

void LifetimeMomentFilter::from_xml(pugi::xml_node node)
{
  this->set_order(std::stoi(get_node_value(node, "order")));
}

void LifetimeMomentFilter::set_order(int order)
{
  if (order < 0) {
    throw std::invalid_argument {"Lifetime moment order must be non-negative."};
  }
  order_ = order;
  n_bins_ = order_ + 1;
}

void LifetimeMomentFilter::get_all_bins(
  const Particle& p, TallyEstimator estimator, FilterMatch& match) const
{
  // Time since the particle was started from its source site. It is reset
  // whenever a particle is initialized from a site (primary, secondary, or
  // split particle) and advanced with the particle clock along each track.
  double tau = p.lifetime();

  // The weight of bin n is tau^n, built by repeated multiplication so that the
  // zeroth moment has a weight of exactly one
  double weight = 1.0;
  for (int n = 0; n < n_bins_; ++n) {
    match.bins_.push_back(n);
    match.weights_.push_back(weight);
    weight *= tau;
  }
}

void LifetimeMomentFilter::to_statepoint(hid_t filter_group) const
{
  Filter::to_statepoint(filter_group);
  write_dataset(filter_group, "order", order_);
}

std::string LifetimeMomentFilter::text_label(int bin) const
{
  return "Lifetime moment, tau^" + std::to_string(bin);
}

//==============================================================================
// C-API functions
//==============================================================================

extern "C" int openmc_lifetime_moment_filter_get_order(
  int32_t index, int* order)
{
  // Make sure this is a valid index to an allocated filter.
  if (int err = verify_filter(index))
    return err;

  // Get a pointer to the filter and downcast.
  const auto& filt_base = model::tally_filters[index].get();
  auto* filt = dynamic_cast<LifetimeMomentFilter*>(filt_base);

  // Check the filter type.
  if (!filt) {
    set_errmsg("Not a lifetime moment filter.");
    return OPENMC_E_INVALID_TYPE;
  }

  // Output the order.
  *order = filt->order();
  return 0;
}

extern "C" int openmc_lifetime_moment_filter_set_order(int32_t index, int order)
{
  // Make sure this is a valid index to an allocated filter.
  if (int err = verify_filter(index))
    return err;

  // Get a pointer to the filter and downcast.
  const auto& filt_base = model::tally_filters[index].get();
  auto* filt = dynamic_cast<LifetimeMomentFilter*>(filt_base);

  // Check the filter type.
  if (!filt) {
    set_errmsg("Not a lifetime moment filter.");
    return OPENMC_E_INVALID_TYPE;
  }

  // Update the filter.
  try {
    filt->set_order(order);
  } catch (const std::invalid_argument& e) {
    set_errmsg(e.what());
    return OPENMC_E_INVALID_ARGUMENT;
  }
  return 0;
}

} // namespace openmc
