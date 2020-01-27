#ifndef OPENMC_TALLIES_FILTER_CELLFROM_H
#define OPENMC_TALLIES_FILTER_CELLFROM_H

#include <string>

#include "openmc/tallies/filter_cell.h"

namespace openmc {

//==============================================================================
//! Specifies which geometric cells particles exit when crossing a surface.
//==============================================================================

class CellFromFilter : public CellFilter
{
public:
  //----------------------------------------------------------------------------
  // Methods

  std::string type() const override {return "cellfrom";}

  void get_all_bins(const Particle* p, TallyEstimator estimator, FilterMatch& match)
  const override;

  std::string text_label(int bin) const override;
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_CELLFROM_H
