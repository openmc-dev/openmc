#ifndef OPENMC_TALLIES_FILTER_AZIMUTHAL_H
#define OPENMC_TALLIES_FILTER_AZIMUTHAL_H

#include "openmc/vector.h"
#include <string>

#include "openmc/span.h"
#include "openmc/tallies/filter.h"

namespace openmc {

//==============================================================================
//! Bins the incident neutron azimuthal angle (relative to the global xy-plane).
//==============================================================================

class AzimuthalFilter : public Filter {
public:
  //----------------------------------------------------------------------------
  // Constructors, destructors

  ~AzimuthalFilter() = default;

  //----------------------------------------------------------------------------
  // Methods

  std::string type_str() const override { return "azimuthal"; }
  FilterType type() const override { return FilterType::AZIMUTHAL; }

  void from_xml(pugi::xml_node node) override;

  void get_all_bins(const Particle& p, TallyEstimator estimator,
    FilterMatch& match) const override;

  void to_statepoint(hid_t filter_group) const override;

  std::string text_label(int bin) const override;

  //----------------------------------------------------------------------------
  // Accessors

  void set_bins(span<double> bins);

private:
  //----------------------------------------------------------------------------
  // Data members

  vector<double> bins_;
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_AZIMUTHAL_H
