#ifndef OPENMC_TALLIES_FILTER_SURFACE_H
#define OPENMC_TALLIES_FILTER_SURFACE_H

#include <cstdint>
#include <unordered_map>
#include <vector>

#include "openmc/tallies/filter.h"

namespace openmc {

//==============================================================================
//! Specifies which surface particles are crossing
//==============================================================================

class SurfaceFilter : public Filter
{
public:
  ~SurfaceFilter() = default;

  std::string type() const override {return "surface";}

  void from_xml(pugi::xml_node node) override;

  void initialize() override;

  void get_all_bins(Particle* p, int estimator, FilterMatch& match)
  const override;

  void to_statepoint(hid_t filter_group) const override;

  std::string text_label(int bin) const override;

  std::vector<int32_t> surfaces_;
  std::unordered_map<int32_t, int> map_;
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_SURFACE_H
