#ifndef OPENMC_TALLIES_FILTER_MATERIAL_H
#define OPENMC_TALLIES_FILTER_MATERIAL_H

#include <cstdint>
#include <unordered_map>
#include <vector>

#include "openmc/tallies/filter.h"

namespace openmc {

//==============================================================================
//! Specifies which material tally events reside in.
//==============================================================================

class MaterialFilter : public Filter
{
public:
  ~MaterialFilter() = default;

  std::string type() const override {return "material";}

  void from_xml(pugi::xml_node node) override;

  void initialize() override;

  void get_all_bins(const Particle* p, int estimator, FilterMatch& match)
  const override;

  void to_statepoint(hid_t filter_group) const override;

  std::string text_label(int bin) const override;

  //! The indices of the materials binned by this filter.
  std::vector<int32_t> materials_;

  //! A map from material indices to filter bin indices.
  std::unordered_map<int32_t, int> map_;
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_MATERIAL_H
