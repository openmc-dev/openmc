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
  ~LegendreFilter() = default;

  std::string type() const override {return "legendre";}

  void from_xml(pugi::xml_node node) override;

  void get_all_bins(const Particle* p, int estimator, FilterMatch& match)
  const override;

  void to_statepoint(hid_t filter_group) const override;

  std::string text_label(int bin) const override;

  int order_;
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_LEGENDRE_H
