#ifndef OPENMC_TALLIES_FILTER_FISSION_ENERGY_H
#define OPENMC_TALLIES_FILTER_FISSION_ENERGY_H

#include "openmc/tallies/filter.h"

namespace openmc {

//==============================================================================
//! Multiplies tally scores by an arbitrary function of incident energy
//! described by a piecewise linear-linear interpolation.
//==============================================================================

class FissionYieldsFilter : public Filter {
public:
  //----------------------------------------------------------------------------
  // Constructors, destructors

  ~FissionYieldsFilter() = default;

  //----------------------------------------------------------------------------
  // Methods

  std::string type_str() const override { return "fissionyields"; }
  FilterType type() const override { return FilterType::FISSION_YIELDS; }

  void from_xml(pugi::xml_node node) override;

  void get_all_bins(const Particle& p, TallyEstimator estimator,
    FilterMatch& match) const override;

  void to_statepoint(hid_t filter_group) const override;

  std::string text_label(int bin) const override;

  //----------------------------------------------------------------------------
  // Accessors

  const vector<std::string>& bins() const { return bins_; }

private:
  //----------------------------------------------------------------------------
  // Data members

  vector<std::string> bins_;
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_FISSION_ENERGY_H
