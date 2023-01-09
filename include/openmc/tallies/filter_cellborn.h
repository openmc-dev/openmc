#ifndef OPENMC_TALLIES_FILTER_CELLBORN_H
#define OPENMC_TALLIES_FILTER_CELLBORN_H

#include <string>

#include "openmc/tallies/filter_cell.h"

namespace openmc {

//==============================================================================
//! Specifies which cell the particle was born in.
//==============================================================================

class CellBornFilter : public CellFilter {
public:
  //----------------------------------------------------------------------------
  // Methods

  std::string type_str() const override { return "cellborn"; }
  FilterType type() const override { return FilterType::CELLBORN; }

  void get_all_bins(const Particle& p, TallyEstimator estimator,
    FilterMatch& match) const override;

  std::string text_label(int bin) const override;
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_CELLBORN_H
