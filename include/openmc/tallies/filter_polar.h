#ifndef OPENMC_TALLIES_FILTER_POLAR_H
#define OPENMC_TALLIES_FILTER_POLAR_H

#include <cmath>

#include "openmc/span.h"
#include "openmc/tallies/filter.h"
#include "openmc/vector.h"

namespace openmc {

//==============================================================================
//! Bins the incident neutron polar angle (relative to the global z-axis).
//==============================================================================

class PolarFilter : public Filter {
public:
  //----------------------------------------------------------------------------
  // Constructors, destructors

  ~PolarFilter() = default;

  //----------------------------------------------------------------------------
  // Methods

  std::string type_str() const override { return "polar"; }
  FilterType type() const override { return FilterType::POLAR; }

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
#endif // OPENMC_TALLIES_FILTER_POLAR_H
