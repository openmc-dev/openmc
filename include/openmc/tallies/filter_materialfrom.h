#ifndef OPENMC_TALLIES_FILTER_MATERIALFROM_H
#define OPENMC_TALLIES_FILTER_MATERIALFROM_H

#include <string>

#include "openmc/tallies/filter_material.h"

namespace openmc {

//==============================================================================
//! Specifies which material particles exit when crossing a surface.
//==============================================================================

class MaterialFromFilter : public MaterialFilter {
public:
  //----------------------------------------------------------------------------
  // Methods

  std::string type_str() const override { return "materialfrom"; }
  FilterType type() const override { return FilterType::MATERIALFROM; }

  void get_all_bins(const Particle& p, TallyEstimator estimator,
    FilterMatch& match) const override;

  std::string text_label(int bin) const override;
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_MATERIALFROM_H
