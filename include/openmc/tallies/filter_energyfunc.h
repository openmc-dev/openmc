#ifndef OPENMC_TALLIES_FILTER_ENERGYFUNC_H
#define OPENMC_TALLIES_FILTER_ENERGYFUNC_H

#include <vector>

#include "openmc/tallies/filter.h"

namespace openmc {

//==============================================================================
//! Multiplies tally scores by an arbitrary function of incident energy
//! described by a piecewise linear-linear interpolation.
//==============================================================================

class EnergyFunctionFilter : public Filter
{
public:
  EnergyFunctionFilter()
    : Filter {}
  {
    n_bins_ = 1;
  }

  ~EnergyFunctionFilter() = default;

  std::string type() const override {return "energyfunction";}

  void from_xml(pugi::xml_node node) override;

  void get_all_bins(Particle* p, int estimator, FilterMatch& match)
  const override;

  void to_statepoint(hid_t filter_group) const override;

  std::string text_label(int bin) const override;

  std::vector<double> energy_;
  std::vector<double> y_;
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_ENERGYFUNC_H
