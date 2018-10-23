#ifndef OPENMC_TALLIES_FILTER_AZIMUTHAL_H
#define OPENMC_TALLIES_FILTER_AZIMUTHAL_H

#include <string>
#include <vector>

#include "openmc/tallies/filter.h"

namespace openmc {

//==============================================================================
//! Bins the incident neutron azimuthal angle (relative to the global xy-plane).
//==============================================================================

class AzimuthalFilter : public Filter
{
public:
  ~AzimuthalFilter() = default;

  std::string type() const override {return "azimuthal";}

  void from_xml(pugi::xml_node node) override;

  void get_all_bins(const Particle* p, int estimator, FilterMatch& match)
  const override;

  void to_statepoint(hid_t filter_group) const override;

  std::string text_label(int bin) const override;

  std::vector<double> bins_;
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_AZIMUTHAL_H
