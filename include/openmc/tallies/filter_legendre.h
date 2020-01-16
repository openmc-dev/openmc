#ifndef OPENMC_TALLIES_FILTER_LEGENDRE_H
#define OPENMC_TALLIES_FILTER_LEGENDRE_H

#include <string>

#include "openmc/tallies/filter.h"

namespace openmc {

//==============================================================================
//! Gives Legendre moments of the change in scattering angle
//==============================================================================

class LegendreFilter : public Filter
{
public:
  //----------------------------------------------------------------------------
  // Constructors, destructors

  ~LegendreFilter() = default;

  //----------------------------------------------------------------------------
  // Methods

  std::string type() const override {return "legendre";}

  void from_xml(pugi::xml_node node) override;

  void get_all_bins(const Particle* p, TallyEstimator estimator, FilterMatch& match)
  const override;

  void to_statepoint(hid_t filter_group) const override;

  std::string text_label(int bin) const override;

  //----------------------------------------------------------------------------
  // Accessors

  int order() const { return order_; }

  void set_order(int order);

private:
  //----------------------------------------------------------------------------
  // Data members

  int order_;
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_LEGENDRE_H
