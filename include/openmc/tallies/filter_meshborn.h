#ifndef OPENMC_TALLIES_FILTER_MESHBORN_H
#define OPENMC_TALLIES_FILTER_MESHBORN_H

#include <cstdint>

#include "openmc/position.h"
#include "openmc/tallies/filter_mesh.h"

namespace openmc {

class MeshBornFilter : public MeshFilter {
public:
  //----------------------------------------------------------------------------
  // Methods

  std::string type_str() const override { return "meshborn"; }
  FilterType type() const override { return FilterType::MESHBORN; }

  void get_all_bins(const Particle& p, TallyEstimator estimator,
    FilterMatch& match) const override;

  std::string text_label(int bin) const override;
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_MESHBORN_H
